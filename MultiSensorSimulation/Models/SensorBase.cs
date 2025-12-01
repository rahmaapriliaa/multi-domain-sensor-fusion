using System;
using System.Numerics;

namespace MultiSensorSimulation.Models
{
    public abstract class SensorBase
    {
        public string? Name { get; set; }
        public string? Unit { get; set; }

        // Gain dan Tau untuk analisis sistem
        public double Gain { get; set; } = 1.0;
        public double Tau { get; set; } = 1.0;

        // Random generator
        protected Random random = new Random();

        // ============================================================
        // UTILITY FUNCTIONS
        // ============================================================

        /// <summary>
        /// Gaussian noise generator menggunakan Box-Muller transform
        /// </summary>
        protected double GaussianNoise(double stdDev)
        {
            double u1 = 1.0 - random.NextDouble();
            double u2 = 1.0 - random.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            return randStdNormal * stdDev;
        }

        /// <summary>
        /// Sinusoidal signal generator
        /// </summary>
        protected double Sinus(double amplitude, double freq, double t)
        {
            return amplitude * Math.Sin(2 * Math.PI * freq * t);
        }

        /// <summary>
        /// ADC quantization simulation
        /// </summary>
        protected double Quantize(double value, double min, double max, int bits)
        {
            double steps = Math.Pow(2, bits);
            double stepSize = (max - min) / steps;

            if (value < min) value = min;
            if (value > max) value = max;

            return Math.Round((value - min) / stepSize) * stepSize + min;
        }

        /// <summary>
        /// Random-walk drift untuk sensor drift simulation
        /// </summary>
        protected double ApplyDrift(double currentDrift, double driftRate, double dt)
        {
            double noise = (random.NextDouble() - 0.5) * 2.0;
            return currentDrift + (noise * driftRate * dt);
        }

        /// <summary>
        /// First-order low-pass filter
        /// </summary>
        protected double LowPassFilter(double input, double prevOutput, double dt)
        {
            double alpha = dt / (Tau + dt);
            return prevOutput + alpha * (input - prevOutput);
        }

        // ============================================================
        // ABSTRACT METHODS - HARUS DIIMPLEMENTASI SEMUA SENSOR
        // ============================================================

        /// <summary>
        /// Generate sinyal sensor dalam time domain
        /// </summary>
        public abstract double[] GenerateSignal(double gain, double tau, double duration);

        /// <summary>
        /// Frequency response untuk Bode plot
        /// </summary>
        public abstract Complex[] GetFrequencyResponse(int points = 1000);

        /// <summary>
        /// Poles dan zeros di S-domain (Laplace)
        /// </summary>
        public abstract Complex[] GetLaplacePolesZeros();

        /// <summary>
        /// Poles dan zeros di Z-domain (discrete)
        /// </summary>
        public abstract Complex[] GetZDomainPolesZeros(double samplingFreq);

        // ============================================================
        // FFT IMPLEMENTATION
        // ============================================================

        public virtual Complex[] CalculateFFT(double[] signal)
        {
            int n = signal.Length;
            int fftSize = 1;
            while (fftSize < n) fftSize *= 2;

            Complex[] fftData = new Complex[fftSize];

            for (int i = 0; i < n; i++)
                fftData[i] = new Complex(signal[i], 0);

            for (int i = n; i < fftSize; i++)
                fftData[i] = Complex.Zero;

            FFT(fftData, false);
            return fftData;
        }

        protected void FFT(Complex[] buffer, bool inverse)
        {
            int n = buffer.Length;
            if (n <= 1) return;

            int bits = (int)Math.Log(n, 2);

            // Bit-reversal permutation
            for (int i = 0; i < n; i++)
            {
                int rev = BitReverse(i, bits);
                if (rev > i)
                {
                    Complex temp = buffer[i];
                    buffer[i] = buffer[rev];
                    buffer[rev] = temp;
                }
            }

            // Cooley-Tukey FFT
            for (int len = 2; len <= n; len *= 2)
            {
                double angle = (inverse ? 2 : -2) * Math.PI / len;
                Complex wlen = new Complex(Math.Cos(angle), Math.Sin(angle));

                for (int i = 0; i < n; i += len)
                {
                    Complex w = Complex.One;
                    for (int j = 0; j < len / 2; j++)
                    {
                        Complex u = buffer[i + j];
                        Complex v = buffer[i + j + len / 2] * w;
                        buffer[i + j] = u + v;
                        buffer[i + j + len / 2] = u - v;
                        w *= wlen;
                    }
                }
            }

            if (inverse)
            {
                for (int i = 0; i < n; i++)
                    buffer[i] /= n;
            }
        }

        private int BitReverse(int n, int bits)
        {
            int reversed = 0;
            for (int i = 0; i < bits; i++)
            {
                reversed = (reversed << 1) | (n & 1);
                n >>= 1;
            }
            return reversed;
        }

        // ============================================================
        // SIGNAL METRICS
        // ============================================================

        public virtual double CalculateRMS(double[] signal)
        {
            double sum = 0;
            foreach (var s in signal) sum += s * s;
            return Math.Sqrt(sum / signal.Length);
        }

        public virtual double CalculatePeakToPeak(double[] signal)
        {
            double min = double.MaxValue, max = double.MinValue;
            foreach (var s in signal)
            {
                if (s < min) min = s;
                if (s > max) max = s;
            }
            return max - min;
        }

        public virtual double CalculateMean(double[] signal)
        {
            double sum = 0;
            foreach (var s in signal) sum += s;
            return sum / signal.Length;
        }

        public virtual double CalculateVariance(double[] signal)
        {
            double mean = CalculateMean(signal);
            double sum = 0;
            foreach (var s in signal)
                sum += Math.Pow(s - mean, 2);
            return sum / signal.Length;
        }

        public virtual double CalculateStdDev(double[] signal)
        {
            return Math.Sqrt(CalculateVariance(signal));
        }
    }
}