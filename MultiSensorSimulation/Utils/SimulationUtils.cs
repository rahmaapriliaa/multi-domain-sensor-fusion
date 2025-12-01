using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;

namespace MultiSensorSimulation.Utils
{
    public static class SimulationUtils
    {
        // ========================================
        // SIGNAL PROCESSING
        // ========================================

        /// <summary>
        /// Normalize data ke range [-1, 1]
        /// </summary>
        public static double[] Normalize(double[] data)
        {
            if (data == null || data.Length == 0)
                return Array.Empty<double>();

            double max = data.Max(Math.Abs);
            return max == 0 ? data : data.Select(x => x / max).ToArray();
        }

        /// <summary>
        /// Normalize ke custom range [min, max]
        /// </summary>
        public static double[] NormalizeToRange(double[] data, double min, double max)
        {
            if (data == null || data.Length == 0)
                return Array.Empty<double>();

            double dataMin = data.Min();
            double dataMax = data.Max();
            double range = dataMax - dataMin;

            if (range == 0)
                return data.Select(x => (min + max) / 2).ToArray();

            return data.Select(x => min + (x - dataMin) / range * (max - min)).ToArray();
        }

        // ========================================
        // FFT OPERATIONS
        // ========================================

        /// <summary>
        /// Compute FFT magnitude dan frequency array
        /// </summary>
        public static (double[] freq, double[] mag) ComputeFFT(double[] signal, double fs)
        {
            if (signal == null || signal.Length == 0)
                return (Array.Empty<double>(), Array.Empty<double>());

            int n = signal.Length;

            // Convert double to Complex32 for MathNet.Numerics
            Complex32[] data = Array.ConvertAll(signal, v => new Complex32((float)v, 0));

            // Forward FFT
            Fourier.Forward(data, FourierOptions.Matlab);

            // Take only positive frequencies (first half)
            int halfN = n / 2;
            double[] mag = data.Take(halfN)
                              .Select(c => (double)c.Magnitude)
                              .ToArray();

            double[] freq = Enumerable.Range(0, halfN)
                                      .Select(i => i * fs / n)
                                      .ToArray();

            return (freq, mag);
        }

        /// <summary>
        /// Compute power spectral density (PSD)
        /// </summary>
        public static (double[] freq, double[] psd) ComputePSD(double[] signal, double fs)
        {
            var (freq, mag) = ComputeFFT(signal, fs);

            // PSD = |FFT|² / N
            double[] psd = mag.Select(m => m * m / signal.Length).ToArray();

            return (freq, psd);
        }

        // ========================================
        // FILTERING
        // ========================================

        /// <summary>
        /// Simple moving average filter
        /// </summary>
        public static double[] MovingAverage(double[] data, int windowSize)
        {
            if (data == null || data.Length == 0 || windowSize <= 0)
                return data;

            double[] filtered = new double[data.Length];

            for (int i = 0; i < data.Length; i++)
            {
                int start = Math.Max(0, i - windowSize / 2);
                int end = Math.Min(data.Length, i + windowSize / 2 + 1);

                double sum = 0;
                for (int j = start; j < end; j++)
                    sum += data[j];

                filtered[i] = sum / (end - start);
            }

            return filtered;
        }

        /// <summary>
        /// First-order low-pass filter (exponential smoothing)
        /// </summary>
        public static double[] LowPassFilter(double[] data, double alpha)
        {
            if (data == null || data.Length == 0)
                return data;

            double[] filtered = new double[data.Length];
            filtered[0] = data[0];

            for (int i = 1; i < data.Length; i++)
            {
                filtered[i] = alpha * data[i] + (1 - alpha) * filtered[i - 1];
            }

            return filtered;
        }

        // ========================================
        // TIME VECTOR GENERATION
        // ========================================

        /// <summary>
        /// Generate time vector untuk plotting
        /// </summary>
        public static double[] GenerateTime(double duration, double dt)
        {
            int n = (int)(duration / dt);
            return Enumerable.Range(0, n).Select(i => i * dt).ToArray();
        }

        /// <summary>
        /// Generate time vector dari sampling rate
        /// </summary>
        public static double[] GenerateTimeFromFs(double duration, double fs)
        {
            return GenerateTime(duration, 1.0 / fs);
        }

        // ========================================
        // SIGNAL COMBINATION
        // ========================================

        /// <summary>
        /// Sum multiple sensor arrays (weighted average)
        /// </summary>
        public static double[] WeightedSum(double[][] arrays, double[] weights)
        {
            if (arrays == null || arrays.Length == 0)
                return Array.Empty<double>();

            if (weights == null || weights.Length != arrays.Length)
            {
                // Default equal weights
                weights = Enumerable.Repeat(1.0 / arrays.Length, arrays.Length).ToArray();
            }

            int n = arrays.Min(a => a.Length);
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < arrays.Length; j++)
                    sum += arrays[j][i] * weights[j];
                result[i] = sum;
            }

            return result;
        }

        /// <summary>
        /// Simple unweighted sum
        /// </summary>
        public static double[] Sum(params double[][] arrays)
        {
            if (arrays == null || arrays.Length == 0)
                return Array.Empty<double>();

            int n = arrays.Min(a => a.Length);
            double[] sum = new double[n];

            foreach (var a in arrays)
                for (int i = 0; i < n; i++)
                    sum[i] += a[i];

            return sum;
        }

        // ========================================
        // STATISTICS
        // ========================================

        /// <summary>
        /// Calculate correlation coefficient antara dua sinyal
        /// </summary>
        public static double CalculateCorrelation(double[] x, double[] y)
        {
            if (x.Length != y.Length)
                throw new ArgumentException("Arrays must have same length");

            int n = x.Length;
            double meanX = x.Average();
            double meanY = y.Average();

            double numerator = 0;
            double denomX = 0;
            double denomY = 0;

            for (int i = 0; i < n; i++)
            {
                double dx = x[i] - meanX;
                double dy = y[i] - meanY;

                numerator += dx * dy;
                denomX += dx * dx;
                denomY += dy * dy;
            }

            return numerator / Math.Sqrt(denomX * denomY);
        }

        /// <summary>
        /// Calculate SNR (Signal-to-Noise Ratio) in dB
        /// </summary>
        public static double CalculateSNR(double[] signal, double[] noise)
        {
            double signalPower = signal.Select(x => x * x).Average();
            double noisePower = noise.Select(x => x * x).Average();

            if (noisePower == 0)
                return double.PositiveInfinity;

            return 10 * Math.Log10(signalPower / noisePower);
        }

        // ========================================
        // WINDOW FUNCTIONS
        // ========================================

        /// <summary>
        /// Hanning window
        /// </summary>
        public static double[] HanningWindow(int length)
        {
            double[] window = new double[length];
            for (int i = 0; i < length; i++)
            {
                window[i] = 0.5 * (1 - Math.Cos(2 * Math.PI * i / (length - 1)));
            }
            return window;
        }

        /// <summary>
        /// Hamming window
        /// </summary>
        public static double[] HammingWindow(int length)
        {
            double[] window = new double[length];
            for (int i = 0; i < length; i++)
            {
                window[i] = 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (length - 1));
            }
            return window;
        }

        /// <summary>
        /// Apply window to signal
        /// </summary>
        public static double[] ApplyWindow(double[] signal, double[] window)
        {
            if (signal.Length != window.Length)
                throw new ArgumentException("Signal and window must have same length");

            return signal.Zip(window, (s, w) => s * w).ToArray();
        }

        // ========================================
        // RESAMPLING
        // ========================================

        /// <summary>
        /// Simple linear interpolation resampling
        /// </summary>
        public static double[] Resample(double[] data, int newLength)
        {
            if (data == null || data.Length == 0 || newLength <= 0)
                return Array.Empty<double>();

            double[] resampled = new double[newLength];
            double ratio = (double)(data.Length - 1) / (newLength - 1);

            for (int i = 0; i < newLength; i++)
            {
                double pos = i * ratio;
                int idx = (int)pos;
                double frac = pos - idx;

                if (idx + 1 < data.Length)
                    resampled[i] = data[idx] * (1 - frac) + data[idx + 1] * frac;
                else
                    resampled[i] = data[idx];
            }

            return resampled;
        }

        // ========================================
        // SIGNAL GENERATION
        // ========================================

        /// <summary>
        /// Generate white noise
        /// </summary>
        public static double[] GenerateWhiteNoise(int length, double amplitude, Random rand = null)
        {
            rand = rand ?? new Random();
            return Enumerable.Range(0, length)
                             .Select(_ => (rand.NextDouble() - 0.5) * 2.0 * amplitude)
                             .ToArray();
        }

        /// <summary>
        /// Generate sine wave
        /// </summary>
        public static double[] GenerateSine(double frequency, double amplitude, double duration, double fs)
        {
            int n = (int)(duration * fs);
            double[] signal = new double[n];

            for (int i = 0; i < n; i++)
            {
                double t = i / fs;
                signal[i] = amplitude * Math.Sin(2 * Math.PI * frequency * t);
            }

            return signal;
        }

        // ========================================
        // VALIDATION
        // ========================================

        /// <summary>
        /// Check if signal contains NaN or Infinity
        /// </summary>
        public static bool IsValid(double[] signal)
        {
            return signal != null &&
                   signal.All(x => !double.IsNaN(x) && !double.IsInfinity(x));
        }

        /// <summary>
        /// Replace NaN/Inf with zero
        /// </summary>
        public static double[] CleanSignal(double[] signal)
        {
            return signal.Select(x =>
                double.IsNaN(x) || double.IsInfinity(x) ? 0.0 : x
            ).ToArray();
        }
    }
}