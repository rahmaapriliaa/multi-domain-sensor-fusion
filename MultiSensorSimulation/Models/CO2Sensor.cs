using System;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public class CO2Sensor : SensorBase
    {
        public double Baseline { get; set; } = 450.0;          // ppm
        public double HumanRise { get; set; } = 300.0;         // tambahan ppm jika ada manusia
        public double TimeConstantResp { get; set; } = 10.0;   // detik
        public double BreathAmplitude { get; set; } = 20.0;    // napas ±20 ppm
        public double BreathFreq { get; set; } = 0.25;         // 0.25 Hz (15 bpm)
        public double NoiseAmplitude { get; set; } = 5.0;      // ±5 ppm
        public double DriftRate { get; set; } = 0.05;          // drift kecil
        public double SamplingRate { get; set; } = 10.0;       // 10 Hz (CO₂ sensor lambat)

        private double drift = 0;

        public override double[] GenerateSignal(double gain, double tau, double duration)
        {
            double plotFs = 100.0;
            int samples = (int)(duration * plotFs);
            double dt = 1.0 / plotFs;
            double[] signal = new double[samples];

            double sensorUpdateInterval = 1.0 / SamplingRate;
            double timeSinceUpdate = 0;
            double lastValue = Baseline;

            for (int i = 0; i < samples; i++)
            {
                double t = i * dt;
                timeSinceUpdate += dt;

                // 1. CO₂ increase due to human (exponential rise)
                double humanCO2 = Baseline + HumanRise * (1 - Math.Exp(-t / TimeConstantResp));

                // 2. Breathing oscillation
                double breathing = BreathAmplitude * Math.Sin(2 * Math.PI * BreathFreq * t);

                // 3. Drift (random walk)
                drift = ApplyDrift(drift, DriftRate, dt);

                // Continuous raw signal
                double raw = humanCO2 + breathing + drift;

                // Quantize + noise only at sampling moments
                if (timeSinceUpdate >= sensorUpdateInterval)
                {
                    double noise = (random.NextDouble() - 0.5) * NoiseAmplitude * 2.0;
                    double sensed = raw + noise;

                    // quantization 1 ppm typical
                    lastValue = Quantize(sensed, 0, 5000, 12);

                    timeSinceUpdate = 0;
                }

                // Sample & hold for plotting rate
                signal[i] = lastValue;
            }

            return signal;
        }

        public override Complex[] GetFrequencyResponse(int points = 1000)
        {
            Complex[] response = new Complex[points];
            double maxFreq = 2.0;

            for (int i = 0; i < points; i++)
            {
                double freq = (i * maxFreq) / points;
                double ω = 2 * Math.PI * freq;
                Complex s = new Complex(0, ω);

                // Low-pass response of CO₂ sensor
                response[i] = 1.0 / (TimeConstantResp * s + 1);
            }

            return response;
        }

        public override Complex[] GetLaplacePolesZeros()
        {
            return new Complex[] { new Complex(-1.0 / TimeConstantResp, 0) };
        }

        public override Complex[] GetZDomainPolesZeros(double samplingFreq)
        {
            double dt = 1.0 / samplingFreq;
            double pole = Math.Exp(-dt / TimeConstantResp);

            return new Complex[]
            {
                new Complex(pole, 0),
                new Complex(0, 0)
            };
        }
    }
}
