using System;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public class MicSensor : SensorBase
    {
        // Properties yang sudah ada
        public double SignalAmplitude { get; set; } = 0.5;
        public double SignalFrequency { get; set; } = 200;
        public double HumAmplitude { get; set; } = 0.05;
        public double HumFrequency { get; set; } = 50;
        public double DetectionThreshold { get; set; } = 0.1;

        // Constructor untuk set Name dan Unit
        public MicSensor()
        {
            Name = "MEMS Microphone";
            Unit = "V";
        }

        public override double[] GenerateSignal(double gain, double tau, double duration)
        {
            int samples = (int)(duration * 100.0);
            double dt = 1.0 / 100.0;
            double[] signal = new double[samples];

            for (int i = 0; i < samples; i++)
            {
                double t = i * dt;
                double voice = SignalAmplitude * gain * Math.Sin(2 * Math.PI * SignalFrequency * t);
                double hum = HumAmplitude * Math.Sin(2 * Math.PI * HumFrequency * t);
                signal[i] = voice + hum;
            }
            return signal;
        }

        public override Complex[] GetFrequencyResponse(int points = 1000)
        {
            Complex[] response = new Complex[points];
            double maxFreq = 10000.0;

            for (int i = 0; i < points; i++)
            {
                double freq = (i * maxFreq) / points;
                response[i] = new Complex(1.0, 0);
            }

            return response;
        }

        public override Complex[] GetLaplacePolesZeros()
        {
            return new Complex[] { new Complex(-1000.0, 0) };
        }

        public override Complex[] GetZDomainPolesZeros(double samplingFreq)
        {
            double dt = 1.0 / samplingFreq;
            double pole = Math.Exp(-1000.0 * dt);

            return new Complex[]
            {
                new Complex(pole, 0),
                new Complex(0, 0)
            };
        }
    }
}