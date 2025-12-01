using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using ScottPlot;

namespace MultiSensorSimulation.Utils
{
    public static class SystemAnalyzer
    {
        // ========================================
        // BASIC TRANSFER FUNCTION MODELS
        // ========================================

        /// <summary>
        /// First-order transfer function: G(s) = K / (τs + 1)
        /// </summary>
        public static (double[] num, double[] den) FirstOrder(double K, double tau)
        {
            return (new double[] { K }, new double[] { tau, 1.0 });
        }

        /// <summary>
        /// Second-order transfer function: G(s) = K*ωn² / (s² + 2ζωn*s + ωn²)
        /// </summary>
        public static (double[] num, double[] den) SecondOrder(double K, double wn, double zeta)
        {
            double[] num = new double[] { K * wn * wn };
            double[] den = new double[] { 1.0, 2.0 * zeta * wn, wn * wn };
            return (num, den);
        }

        // ========================================
        // EXTRACT TRANSFER FUNCTION FROM SENSORS
        // ========================================

        /// <summary>
        /// Extract transfer function parameters from sensor instance
        /// </summary>
        public static (double[] num, double[] den) GetSensorTransferFunction(Models.SensorBase sensor)
        {
            // Default first-order system
            double gain = 1.0;
            double tau = 0.1;

            // MLX90640 Thermal Sensor
            if (sensor is Models.MLX90640Sensor thermal)
            {
                gain = thermal.Emissivity; // Temperature response gain
                tau = thermal.TimeConstant; // Thermal time constant
                return FirstOrder(gain, tau);
            }

            // CO2 Sensor
            else if (sensor is Models.CO2Sensor co2)
            {
                gain = co2.HumanRise / co2.Baseline; // Normalized response
                tau = co2.TimeConstantResp; // Respiratory response time
                return FirstOrder(gain, tau);
            }

            // Radar Sensor (Second-order bandpass)
            else if (sensor is Models.RadarSensor radar)
            {
                double dopplerFreq = (2 * radar.TargetVelocity * radar.CarrierFrequency) / 3.0e8;
                double wn = 2 * Math.PI * dopplerFreq;
                double zeta = 0.7; // Typical damping for radar

                // Bandpass approximation: s / (s² + 2ζωn*s + ωn²)
                double[] num = new double[] { radar.DopplerGain, 0 };
                double[] den = new double[] { 1.0, 2.0 * zeta * wn, wn * wn };
                return (num, den);
            }

            // IMU Sensor (High-pass + Low-pass cascade)
            else if (sensor is Models.IMUSensor imu)
            {
                // Simplified low-pass model for IMU
                double cutoff = 40.0; // Hz
                tau = 1.0 / (2 * Math.PI * cutoff);
                gain = imu.AccelSensitivity;
                return FirstOrder(gain, tau);
            }

            // Microphone (Flat response with high-pass)
            else if (sensor is Models.MicSensor mic)
            {
                gain = mic.SignalAmplitude;
                tau = 0.001; // Very fast response (1ms)
                return FirstOrder(gain, tau);
            }

            // Default fallback
            return FirstOrder(gain, tau);
        }

        // ========================================
        // COMBINE MULTIPLE TRANSFER FUNCTIONS
        // ========================================

        /// <summary>
        /// Combine multiple transfer functions (parallel addition)
        /// Used for sensor fusion analysis
        /// </summary>
        public static (double[] num, double[] den) Combine(List<(double[] num, double[] den)> list)
        {
            if (list == null || list.Count == 0)
                return (new double[] { 0.0 }, new double[] { 1.0 });

            // Find common denominator (LCM of all denominators)
            double[] commonDen = new double[] { 1.0 };
            foreach (var (_, den) in list)
            {
                commonDen = PolyMultiply(commonDen, den);
            }

            // Sum all numerators scaled by appropriate factors
            double[] totalNum = new double[commonDen.Length];
            foreach (var (num, den) in list)
            {
                double[] scaleFactor = PolyDivide(commonDen, den);
                double[] scaledNum = PolyMultiply(num, scaleFactor);
                totalNum = PolyAdd(totalNum, scaledNum);
            }

            return (Trim(totalNum), Trim(commonDen));
        }

        /// <summary>
        /// Combine sensor transfer functions with weights (for fusion analysis)
        /// </summary>
        public static (double[] num, double[] den) CombineSensors(
            List<Models.SensorBase> sensors,
            List<double> weights = null)
        {
            if (sensors == null || sensors.Count == 0)
                return (new double[] { 0.0 }, new double[] { 1.0 });

            // Default equal weights
            if (weights == null || weights.Count != sensors.Count)
            {
                weights = Enumerable.Repeat(1.0 / sensors.Count, sensors.Count).ToList();
            }

            var tfList = new List<(double[] num, double[] den)>();

            for (int i = 0; i < sensors.Count; i++)
            {
                var (num, den) = GetSensorTransferFunction(sensors[i]);

                // Apply weight to numerator
                double[] weightedNum = num.Select(x => x * weights[i]).ToArray();
                tfList.Add((weightedNum, den));
            }

            return Combine(tfList);
        }

        // ========================================
        // DISCRETIZATION METHODS
        // ========================================

        /// <summary>
        /// Bilinear transform (Tustin's method) for s → z transformation
        /// Most accurate for general systems
        /// </summary>
        public static (double[] numZ, double[] denZ) DiscretizeToBilinear(
            double[] numS, double[] denS, double Ts)
        {
            int numOrder = numS.Length - 1;
            int denOrder = denS.Length - 1;
            int maxOrder = Math.Max(numOrder, denOrder);

            var numZ = new double[maxOrder + 1];
            var denZ = new double[maxOrder + 1];

            double c = 2.0 / Ts; // Bilinear constant

            // Apply bilinear substitution: s = (2/Ts) * (z-1)/(z+1)
            for (int i = 0; i <= numOrder; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    int zPower = i - j;
                    double coeff = numS[numOrder - i] * Binomial(i, j) *
                                  Math.Pow(c, i - j) * Math.Pow(-1, j);
                    if (zPower < numZ.Length)
                        numZ[zPower] += coeff;
                }
            }

            for (int i = 0; i <= denOrder; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    int zPower = i - j;
                    double coeff = denS[denOrder - i] * Binomial(i, j) *
                                  Math.Pow(c, i - j) * Math.Pow(-1, j);
                    if (zPower < denZ.Length)
                        denZ[zPower] += coeff;
                }
            }

            // Normalize by leading coefficient
            double lead = denZ[0];
            if (Math.Abs(lead) > 1e-12)
            {
                for (int i = 0; i < denZ.Length; i++) denZ[i] /= lead;
                for (int i = 0; i < numZ.Length; i++) numZ[i] /= lead;
            }

            return (Trim(numZ), Trim(denZ));
        }

        /// <summary>
        /// Matched pole-zero method: z = exp(s*Ts)
        /// Good for systems with known poles/zeros
        /// </summary>
        public static (double[] numZ, double[] denZ) DiscretizeToZ(
            double[] numS, double[] denS, double Ts)
        {
            var polesS = FindRoots(denS);
            var zerosS = FindRoots(numS);

            // Map s-plane to z-plane: z = exp(s*Ts)
            var polesZ = polesS.Select(p => Complex.Exp(p * Ts)).ToArray();
            var zerosZ = zerosS.Select(z => Complex.Exp(z * Ts)).ToArray();

            // Reconstruct polynomials from roots
            var numZc = PolyFromRoots(zerosZ);
            var denZc = PolyFromRoots(polesZ);

            // Extract real parts (should be real for real systems)
            double[] numZ = numZc.Select(c => c.Real).ToArray();
            double[] denZ = denZc.Select(c => c.Real).ToArray();

            // Normalize
            double lead = denZ[0];
            if (Math.Abs(lead) > 1e-12)
            {
                for (int i = 0; i < denZ.Length; i++) denZ[i] /= lead;
                for (int i = 0; i < numZ.Length; i++) numZ[i] /= lead;
            }

            return (Trim(numZ), Trim(denZ));
        }

        // ========================================
        // PLOTTING FUNCTIONS
        // ========================================

        /// <summary>
        /// Plot discrete frequency response (Bode magnitude)
        /// </summary>
        public static void PlotFrequencyResponse(
            Plot plot, double[] num, double[] den, double fs, string title)
        {
            plot.Clear();

            int numPoints = 300;
            double fMin = 0.01;
            double fMax = fs / 2;

            // Logarithmic frequency spacing
            double[] frequencies = Enumerable.Range(0, numPoints)
                .Select(i => fMin * Math.Pow(fMax / fMin, (double)i / (numPoints - 1)))
                .ToArray();

            double[] magnitudeDB = new double[numPoints];

            for (int i = 0; i < numPoints; i++)
            {
                double omega = 2 * Math.PI * frequencies[i] / fs;
                Complex z = Complex.Exp(new Complex(0, omega));

                Complex numer = EvaluatePoly(num, z);
                Complex denom = EvaluatePoly(den, z);
                Complex H = numer / (denom + 1e-12); // Prevent division by zero

                magnitudeDB[i] = 20 * Math.Log10(H.Magnitude + 1e-12);
            }

            var magScatter = plot.Add.Scatter(frequencies, magnitudeDB);
            magScatter.Color = ScottPlot.Colors.Cyan;
            magScatter.LineWidth = 2;

            plot.Title(title);
            plot.Axes.Bottom.Label.Text = "Frequency (Hz)";
            plot.Axes.Left.Label.Text = "Magnitude (dB)";
            plot.Axes.SetLimitsX(fMin, fMax);
        }

        /// <summary>
        /// Plot continuous transfer function (s-domain Bode)
        /// </summary>
        public static void PlotTransferFunction(
            Plot plot, double[] num, double[] den, string title)
        {
            plot.Clear();

            int numPoints = 500;
            double wMin = 0.01;
            double wMax = 1000.0;

            double[] w = Enumerable.Range(0, numPoints)
                .Select(i => wMin * Math.Pow(wMax / wMin, (double)i / numPoints))
                .ToArray();

            double[] mag = new double[numPoints];

            for (int i = 0; i < numPoints; i++)
            {
                var jw = new Complex(0, w[i]);

                Complex nval = EvaluatePoly(num, jw);
                Complex dval = EvaluatePoly(den, jw);
                Complex H = nval / (dval + 1e-12);

                mag[i] = 20 * Math.Log10(H.Magnitude + 1e-12);
            }

            var scatter = plot.Add.Scatter(w, mag);
            scatter.Color = ScottPlot.Colors.White;
            scatter.LineWidth = 2;

            plot.Title(title);
            plot.Axes.Bottom.Label.Text = "ω (rad/s)";
            plot.Axes.Left.Label.Text = "Magnitude (dB)";
            plot.Axes.SetLimitsX(wMin, wMax);
        }

        // ========================================
        // POLYNOMIAL OPERATIONS
        // ========================================

        private static Complex EvaluatePoly(double[] poly, Complex x)
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < poly.Length; i++)
            {
                result += poly[i] * Complex.Pow(x, poly.Length - 1 - i);
            }
            return result;
        }

        private static double Binomial(int n, int k)
        {
            if (k > n || k < 0) return 0;
            if (k == 0 || k == n) return 1;

            double result = 1;
            for (int i = 1; i <= k; i++)
            {
                result *= (n - k + i);
                result /= i;
            }
            return result;
        }

        private static double[] PolyMultiply(double[] a, double[] b)
        {
            double[] r = new double[a.Length + b.Length - 1];
            for (int i = 0; i < a.Length; i++)
                for (int j = 0; j < b.Length; j++)
                    r[i + j] += a[i] * b[j];
            return r;
        }

        private static double[] PolyAdd(double[] a, double[] b)
        {
            int n = Math.Max(a.Length, b.Length);
            double[] A = new double[n];
            double[] B = new double[n];

            Array.Copy(a, 0, A, n - a.Length, a.Length);
            Array.Copy(b, 0, B, n - b.Length, b.Length);

            double[] R = new double[n];
            for (int i = 0; i < n; i++)
                R[i] = A[i] + B[i];

            return R;
        }

        private static double[] PolyDivide(double[] dividend, double[] divisor)
        {
            dividend = (double[])dividend.Clone();
            int n = dividend.Length;
            int m = divisor.Length;

            if (m > n) return new double[] { 0.0 };

            double[] q = new double[n - m + 1];

            for (int k = 0; k <= n - m; k++)
            {
                if (Math.Abs(divisor[0]) < 1e-12)
                    return new double[] { 0.0 };

                double coeff = dividend[k] / divisor[0];
                q[k] = coeff;

                for (int j = 0; j < m; j++)
                    dividend[k + j] -= coeff * divisor[j];
            }

            return Trim(q);
        }

        /// <summary>
        /// Find roots of polynomial using eigenvalue method
        /// </summary>
        public static Complex[] FindRoots(double[] poly)
        {
            int n = poly.Length - 1;
            if (n < 1) return Array.Empty<Complex>();

            double lead = poly[0];
            if (Math.Abs(lead) < 1e-12) return Array.Empty<Complex>();

            var c = poly.Select(x => x / lead).ToArray();

            // Companion matrix method
            var M = DenseMatrix.Create(n, n, Complex.Zero);

            for (int i = 0; i < n - 1; i++)
                M[i + 1, i] = Complex.One;

            for (int j = 0; j < n; j++)
                M[0, j] = new Complex(-c[j + 1], 0);

            var evd = M.Evd();
            return evd.EigenValues.ToArray();
        }

        private static Complex[] PolyFromRoots(Complex[] roots)
        {
            Complex[] poly = new Complex[] { Complex.One };

            foreach (var r in roots)
            {
                Complex[] next = new Complex[poly.Length + 1];
                for (int i = 0; i < poly.Length; i++)
                {
                    next[i] += poly[i];
                    next[i + 1] -= poly[i] * r;
                }
                poly = next;
            }

            return poly;
        }

        private static double[] Trim(double[] p)
        {
            int idx = 0;
            while (idx < p.Length && Math.Abs(p[idx]) < 1e-14)
                idx++;

            if (idx == p.Length)
                return new double[] { 0.0 };

            double[] r = new double[p.Length - idx];
            Array.Copy(p, idx, r, 0, r.Length);
            return r;
        }
    }
}