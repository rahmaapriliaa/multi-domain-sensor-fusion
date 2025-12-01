using System;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public class MLX90640Sensor : SensorBase
    {
        public double Emissivity { get; set; } = 0.95;
        public double TimeConstant { get; set; } = 2.0;     // respons termal ke target
        public double AmbientTemp { get; set; } = 25.0;      // °C
        public double BodyTemp { get; set; } = 37.0;         // °C (human)
        public double NoiseAmplitude { get; set; } = 0.05;   // °C (realistic)
        public double SamplingRate { get; set; } = 16.0;     // MLX90640 typical FPS

        public override double[] GenerateSignal(double gain, double tau, double duration)
        {
            double plotFs = 100.0;
            int samples = (int)(duration * plotFs);
            double dt = 1.0 / plotFs;
            double[] signal = new double[samples];

            // ΔT human vs ambient
            double deltaT = BodyTemp - AmbientTemp;

            double currentTemp = AmbientTemp;
            double lastSampledVal = AmbientTemp;

            double sensorUpdateInterval = 1.0 / SamplingRate;
            double timeSinceLastUpdate = 0;

            for (int i = 0; i < samples; i++)
            {
                double t = i * dt;
                timeSinceLastUpdate += dt;

                // MODEL THERMAL HUMAN SIGNATURE (benar)
                double targetTemp =
                    AmbientTemp +
                    deltaT * (1 - Math.Exp(-t / TimeConstant));

                // first-order lag (fisik sensor)
                currentTemp += (targetTemp - currentTemp) * dt / TimeConstant;

                // Sample & Hold (sesuai frame rate MLX90640)
                if (timeSinceLastUpdate >= sensorUpdateInterval)
                {
                    double noise = (random.NextDouble() - 0.5) * NoiseAmplitude * 2.0;
                    double raw = currentTemp + noise;

                    // Quantization 0.1°C
                    lastSampledVal = Math.Round(raw * 10) / 10.0;

                    timeSinceLastUpdate = 0;
                }

                signal[i] = lastSampledVal;
            }

            return signal;
        }

        public override Complex[] GetFrequencyResponse(int points = 1000)
        {
            Complex[] response = new Complex[points];
            double maxFreq = 10.0;

            for (int i = 0; i < points; i++)
            {
                double freq = (i * maxFreq) / points;
                double omega = 2 * Math.PI * freq;
                Complex s = new Complex(0, omega);

                response[i] = 1.0 / (TimeConstant * s + 1);
            }

            return response;
        }

        public override Complex[] GetLaplacePolesZeros()
        {
            return new Complex[] { new Complex(-1.0 / TimeConstant, 0) };
        }

        public override Complex[] GetZDomainPolesZeros(double samplingFreq)
        {
            double dt = 1.0 / samplingFreq;
            double pole = Math.Exp(-dt / TimeConstant);

            return new Complex[]
            {
                new Complex(pole, 0),  // pole
                new Complex(0, 0)      // zero
            };
        }
    }
}
