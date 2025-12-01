using System;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public class IMUSensor : SensorBase
    {
        // ==================== SLIDER PARAMETERS =====================
        public double AccelBias { get; set; } = 0.05;
        public double GyroDrift { get; set; } = 0.01;

        public double AccelSensitivity { get; set; } = 1.0;
        public double GyroSensitivity { get; set; } = 1.0;

        public double NoiseAccel { get; set; } = 0.02;
        public double NoiseGyro { get; set; } = 0.01;

        public double SamplingRate { get; set; } = 100.0; // Hz
        // ============================================================

        private double gyroBiasState = 0; // Random walk bias

        public override double[] GenerateSignal(double gain, double tau, double duration)
        {
            int samples = (int)(duration * SamplingRate);
            double dt = 1.0 / SamplingRate;

            double[] accel = new double[samples];
            double[] gyro = new double[samples];

            gyroBiasState = (random.NextDouble() - 0.5) * 0.1; // initial bias

            for (int i = 0; i < samples; i++)
            {
                double t = i * dt;

                // =================== ACCELEROMETER =====================
                // Low frequency breathing / micro movement
                double A_motion = AccelSensitivity * gain *
                                  Math.Sin(2 * Math.PI * 0.5 * t);

                double A_noise = (random.NextDouble() - 0.5) * NoiseAccel;

                accel[i] = A_motion + AccelBias + A_noise;


                // ===================== GYROSCOPE ========================
                // Rotational micro-movement
                double G_motion = GyroSensitivity * gain *
                                  Math.Sin(2 * Math.PI * 1.0 * t);

                // Random-walk drift (slow drift)
                gyroBiasState += GyroDrift * Math.Sqrt(dt) *
                                 (random.NextDouble() - 0.5);

                double G_noise = (random.NextDouble() - 0.5) * NoiseGyro;

                gyro[i] = G_motion + gyroBiasState + G_noise;
            }

            // ===================== RETURN MERGED =======================
            // Format: EVEN = ACCEL, ODD = GYRO
            double[] merged = new double[samples * 2];

            for (int i = 0; i < samples; i++)
            {
                merged[2 * i] = accel[i];
                merged[2 * i + 1] = gyro[i];
            }

            return merged;
        }

        // ============================================================
        // FREQUENCY RESPONSE
        // ============================================================
        public override Complex[] GetFrequencyResponse(int points = 1000)
        {
            Complex[] response = new Complex[points];
            double maxFreq = SamplingRate / 2;

            for (int i = 0; i < points; i++)
            {
                double freq = (i * maxFreq) / points;
                double cutoff = 40.0;

                double mag = 1.0 / Math.Sqrt(1 + Math.Pow(freq / cutoff, 4));
                response[i] = new Complex(mag, 0);
            }
            return response;
        }

        public override Complex[] GetLaplacePolesZeros()
        {
            // FIXED: 4th order Butterworth low-pass filter at 40 Hz
            double cutoff = 40.0;
            double ωc = 2 * Math.PI * cutoff;

            // Butterworth poles at 45° and 135° angles
            double θ1 = Math.PI / 4;      // 45°
            double θ2 = 3 * Math.PI / 4;  // 135°

            Complex pole1 = new Complex(-ωc * Math.Cos(θ1), ωc * Math.Sin(θ1));
            Complex pole2 = new Complex(-ωc * Math.Cos(θ1), -ωc * Math.Sin(θ1));
            Complex pole3 = new Complex(-ωc * Math.Cos(θ2), ωc * Math.Sin(θ2));
            Complex pole4 = new Complex(-ωc * Math.Cos(θ2), -ωc * Math.Sin(θ2));

            return new Complex[] { pole1, pole2, pole3, pole4 };
        }

        public override Complex[] GetZDomainPolesZeros(double samplingFreq)
        {
            double dt = 1.0 / samplingFreq;
            double cutoff = 40.0;
            double ωc = 2 * Math.PI * cutoff;

            double θ1 = Math.PI / 4;
            double θ2 = 3 * Math.PI / 4;

            Complex s1 = new Complex(-ωc * Math.Cos(θ1), ωc * Math.Sin(θ1));
            Complex s2 = new Complex(-ωc * Math.Cos(θ2), ωc * Math.Sin(θ2));

            Complex z1 = Complex.Exp(s1 * dt);
            Complex z2 = Complex.Exp(s2 * dt);

            return new Complex[] { z1, Complex.Conjugate(z1), z2, Complex.Conjugate(z2) };
        }
    }
}