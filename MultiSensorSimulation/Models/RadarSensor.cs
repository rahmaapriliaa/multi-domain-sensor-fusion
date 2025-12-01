using System;
using System.Linq;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public class RadarSensor : SensorBase
    {
        // Sensor-specific parameters
        public double CarrierFrequency { get; set; } = 24.0e9; // Hz (24 GHz K-band)
        public double TargetVelocity { get; set; } = 2.0; // m/s (walking speed)
        public double DopplerGain { get; set; } = 1.0;
        public double ClutterLevel { get; set; } = 0.1; // normalized
        public double NoiseVariance { get; set; } = 0.05;
        public double SamplingRate { get; set; } = 1000.0; // Hz (typical radar PRF)

        private const double SpeedOfLight = 3.0e8; // m/s

        public override double[] GenerateSignal(double gain, double tau, double duration)
        {
            int samples = (int)(duration * SamplingRate);
            double dt = 1.0 / SamplingRate;
            double[] signal = new double[samples];

            // Calculate Doppler frequency shift
            // f_d = (2 * v * f_c) / c
            double dopplerFreq = (2 * TargetVelocity * CarrierFrequency) / SpeedOfLight;

            for (int i = 0; i < samples; i++)
            {
                double t = i * dt;

                // Doppler shifted signal (I/Q baseband)
                double dopplerSignal = Math.Cos(2 * Math.PI * dopplerFreq * t) * gain * DopplerGain;

                // Add range-varying amplitude (simple model)
                double rangeEffect = 1.0 / (1.0 + 0.1 * Math.Abs(Math.Sin(0.5 * t)));

                // Clutter (static objects)
                double clutter = ClutterLevel * Math.Sin(2 * Math.PI * 5 * t) *
                               (random.NextDouble() - 0.5);

                // Thermal noise
                double noise = GenerateGaussianNoise() * Math.Sqrt(NoiseVariance);

                // Multipath interference
                double multipath = 0.2 * Math.Sin(2 * Math.PI * dopplerFreq * t + Math.PI / 4);

                signal[i] = dopplerSignal * rangeEffect + clutter + noise + multipath;
            }

            return signal;
        }

        private double GenerateGaussianNoise()
        {
            double u1 = 1.0 - random.NextDouble();
            double u2 = 1.0 - random.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
        }

        public override Complex[] GetFrequencyResponse(int points = 1000)
        {
            Complex[] response = new Complex[points];
            double maxFreq = SamplingRate / 2;

            for (int i = 0; i < points; i++)
            {
                double freq = (i * maxFreq) / points;

                // Doppler radar response (simplified bandpass filter)
                double dopplerFreq = (2 * TargetVelocity * CarrierFrequency) / SpeedOfLight;
                double bandwidth = 50.0; // Hz

                double response_mag = Math.Exp(-Math.Pow((freq - dopplerFreq) / bandwidth, 2));
                response[i] = new Complex(response_mag, 0);
            }

            return response;
        }

        public override Complex[] GetLaplacePolesZeros()
        {
            // FIXED: Return only DISTINCT poles
            // Simplified second-order bandpass
            double centerFreq = (2 * TargetVelocity * CarrierFrequency) / SpeedOfLight;
            double omega = 2 * Math.PI * centerFreq;
            double zeta = 0.7; // damping ratio

            // Return conjugate pair as separate entries
            Complex pole1 = new Complex(-zeta * omega, omega * Math.Sqrt(1 - zeta * zeta));
            Complex pole2 = new Complex(-zeta * omega, -omega * Math.Sqrt(1 - zeta * zeta));

            return new Complex[] { pole1, pole2 };
        }

        public override Complex[] GetZDomainPolesZeros(double samplingFreq)
        {
            double dt = 1.0 / samplingFreq;
            double centerFreq = (2 * TargetVelocity * CarrierFrequency) / SpeedOfLight;
            double omega = 2 * Math.PI * centerFreq;
            double zeta = 0.7;

            Complex s_pole = new Complex(-zeta * omega, omega * Math.Sqrt(1 - zeta * zeta));
            Complex z_pole = Complex.Exp(s_pole * dt);

            // Return conjugate pair
            return new Complex[] { z_pole, Complex.Conjugate(z_pole) };
        }
    }
}