using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Windows;
using WpfColor = System.Windows.Media.Color;
using WpfColors = System.Windows.Media.Colors;
using System.Windows.Media;
using ScottPlot;
using ScottPlot.WPF;
using MultiSensorSimulation.Models;
using MultiSensorSimulation.Utils;

namespace MultiSensorSimulation
{
    public partial class MainWindow : Window
    {
        private DetectionPipeline detectionPipeline = null!;
        private Dictionary<string, double[]> sensorSignals = new Dictionary<string, double[]>();
        private Dictionary<string, Complex[]> sensorFFTs = new Dictionary<string, Complex[]>();
        private Dictionary<string, double> latestConfidences = new Dictionary<string, double>();
        private bool isInitialized = false;

        public MainWindow()
        {
            InitializeComponent();
            InitializeSystem();
            isInitialized = true;
            RunSimulation();
            LogStatus("✅ System initialized successfully");
        }

        private void InitializeSystem()
        {
            detectionPipeline = new DetectionPipeline
            {
                ConfidenceThreshold = 70.0,
                MinimumPositiveSensors = 3
            };
            InitializeAllPlots();
        }

        private void InitializeAllPlots()
        {
            // Time domain plots
            InitPlot(PlotThermalTime, "Time (s)", "Temperature (°C)", ScottPlot.Colors.OrangeRed);
            InitPlot(PlotCO2Time, "Time (s)", "CO2 (ppm)", ScottPlot.Colors.ForestGreen);
            InitPlot(PlotRadarTime, "Time (s)", "Amplitude (V)", ScottPlot.Colors.DodgerBlue);
            InitPlot(PlotMicTime, "Time (s)", "Amplitude (V)", ScottPlot.Colors.Purple);
            InitPlot(PlotIMUTime, "Time (s)", "Accel (m/s²)", ScottPlot.Colors.Orange);

            // Frequency domain plots
            InitPlot(PlotThermalFFT, "Frequency (Hz)", "Magnitude (dB)", ScottPlot.Colors.OrangeRed);
            InitPlot(PlotCO2FFT, "Frequency (Hz)", "Magnitude (dB)", ScottPlot.Colors.ForestGreen);
            InitPlot(PlotRadarFFT, "Frequency (Hz)", "Magnitude (dB)", ScottPlot.Colors.DodgerBlue);
            InitPlot(PlotMicFFT, "Frequency (Hz)", "Magnitude (dB)", ScottPlot.Colors.Purple);
            InitPlot(PlotIMUFFT, "Frequency (Hz)", "Magnitude (dB)", ScottPlot.Colors.Orange);

            // System plots
            InitPlot(PlotLaplace, "Real", "Imaginary", ScottPlot.Colors.Black);
            InitPlot(PlotZDomain, "Real", "Imaginary", ScottPlot.Colors.Black);

            SetupZDomainUnitCircle();
        }

        private void InitPlot(WpfPlot plot, string xlabel, string ylabel, ScottPlot.Color color)
        {
            plot.Plot.Clear();
            plot.Plot.FigureBackground.Color = ScottPlot.Colors.White;
            plot.Plot.DataBackground.Color = ScottPlot.Colors.White;
            plot.Plot.Axes.Color(ScottPlot.Colors.Black);
            plot.Plot.Grid.MajorLineColor = ScottPlot.Colors.LightGray.WithOpacity(0.3);
            plot.Plot.Axes.Bottom.Label.Text = xlabel;
            plot.Plot.Axes.Bottom.Label.FontSize = 11;
            plot.Plot.Axes.Left.Label.Text = ylabel;
            plot.Plot.Axes.Left.Label.FontSize = 11;
            plot.Refresh();
        }

        private void SetupZDomainUnitCircle()
        {
            PlotZDomain.Plot.Clear();
            var circle = PlotZDomain.Plot.Add.Circle(0, 0, 1);
            circle.LineColor = ScottPlot.Colors.Gray;
            circle.LinePattern = LinePattern.Dashed;
            circle.LineWidth = 2;
            circle.FillColor = ScottPlot.Colors.Transparent;

            PlotZDomain.Plot.Add.VerticalLine(0, color: ScottPlot.Colors.Black, width: 1);
            PlotZDomain.Plot.Add.HorizontalLine(0, color: ScottPlot.Colors.Black, width: 1);

            PlotZDomain.Plot.Axes.SquareUnits();
            PlotZDomain.Plot.Axes.SetLimits(-1.5, 1.5, -1.5, 1.5);
            PlotZDomain.Refresh();
        }

        private void OnParameterChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if (!isInitialized) return;
            try
            {
                UpdateSensorParameters();
                RunSimulation();
            }
            catch (Exception ex)
            {
                LogStatus($"⚠️ Parameter update error: {ex.Message}");
            }
        }

        private void BtnRunSimulation_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                RunSimulation();
                LogStatus("🔄 Manual simulation run completed");
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Simulation Error: {ex.Message}", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void BtnReset_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                ResetAllParameters();
                RunSimulation();
                LogStatus("🔄 All parameters reset to default");
            }
            catch (Exception ex)
            {
                LogStatus($"❌ Reset error: {ex.Message}");
            }
        }

        private void RunSimulation()
        {
            try
            {
                double gain = SliderGain.Value;
                double duration = SliderDuration.Value;

                UpdateSensorParameters();

                sensorSignals = detectionPipeline.RunAllSensors(gain, 0.0, duration);
                CalculateAllFFTs();

                // Update Time Plots
                UpdateTimePlot(PlotThermalTime, sensorSignals["Thermal"], detectionPipeline.GetSensor<MLX90640Sensor>().SamplingRate, "#E74C3C");
                UpdateTimePlot(PlotCO2Time, sensorSignals["CO2"], detectionPipeline.GetSensor<CO2Sensor>().SamplingRate, "#27AE60");
                UpdateTimePlot(PlotRadarTime, sensorSignals["Radar"], detectionPipeline.GetSensor<RadarSensor>().SamplingRate, "#3498DB");
                UpdateTimePlot(PlotMicTime, sensorSignals["Microphone"], 100.0, "#9B59B6");
                UpdateTimePlot(PlotIMUTime, sensorSignals["IMU"], detectionPipeline.GetSensor<IMUSensor>().SamplingRate, "#E67E22");

                // Update FFT Plots
                UpdateFFTPlotDiscrete(PlotThermalFFT, sensorFFTs["Thermal"], detectionPipeline.GetSensor<MLX90640Sensor>().SamplingRate, "#E74C3C");
                UpdateFFTPlotDiscrete(PlotCO2FFT, sensorFFTs["CO2"], detectionPipeline.GetSensor<CO2Sensor>().SamplingRate, "#27AE60");
                UpdateFFTPlotDiscrete(PlotRadarFFT, sensorFFTs["Radar"], detectionPipeline.GetSensor<RadarSensor>().SamplingRate, "#3498DB");
                UpdateFFTPlotDiscrete(PlotMicFFT, sensorFFTs["Microphone"], 100.0, "#9B59B6");
                UpdateFFTPlotDiscrete(PlotIMUFFT, sensorFFTs["IMU"], detectionPipeline.GetSensor<IMUSensor>().SamplingRate, "#E67E22");

                var (detected, confidence, confidences) = detectionPipeline.DetectHuman(sensorSignals);
                latestConfidences = confidences;

                UpdateDetectionStatus(detected, confidence, confidences);
                UpdateLaplacePlot();
                UpdateZDomainPlot();
                UpdateSignalMetrics();

                TxtTimestamp.Text = $"Last Update: {DateTime.Now:HH:mm:ss}";
            }
            catch (Exception ex)
            {
                LogStatus($"❌ Simulation error: {ex.Message}");
            }
        }

        private void UpdateSensorParameters()
        {
            var thermal = detectionPipeline.GetSensor<MLX90640Sensor>();
            var co2 = detectionPipeline.GetSensor<CO2Sensor>();
            var radar = detectionPipeline.GetSensor<RadarSensor>();
            var mic = detectionPipeline.GetSensor<MicSensor>();
            var imu = detectionPipeline.GetSensor<IMUSensor>();

            if (thermal != null)
            {
                thermal.TimeConstant = SliderThermalTau.Value;
                thermal.NoiseAmplitude = SliderThermalNoise.Value * 0.1;
                thermal.SamplingRate = SliderThermalRate.Value;
            }

            if (co2 != null)
            {
                co2.TimeConstantResp = SliderCO2Tau.Value;
                co2.HumanRise = SliderCO2Rise.Value;
                co2.BreathAmplitude = SliderCO2Breath.Value * 0.5;
            }

            if (radar != null)
            {
                radar.TargetVelocity = SliderRadarVel.Value;
                radar.DopplerGain = SliderRadarGain.Value;
                radar.ClutterLevel = SliderRadarClutter.Value * 0.3;
            }

            if (mic != null)
            {
                mic.SignalFrequency = SliderMicFreq.Value;
                mic.SignalAmplitude = SliderMicAmp.Value;
                mic.HumAmplitude = SliderMicHum.Value * 0.2;
            }

            if (imu != null)
            {
                imu.AccelBias = SliderIMUBias.Value * 0.5;
                imu.GyroDrift = SliderIMUDrift.Value * 0.5;
                imu.NoiseAccel = SliderIMUNoise.Value * 0.3;
            }
        }

        private void ResetAllParameters()
        {
            SliderGain.Value = 1.0;
            SliderDuration.Value = 10.0;
            SliderThermalTau.Value = 2.0;
            SliderThermalNoise.Value = 0.05;
            SliderThermalRate.Value = 16.0;
            SliderCO2Tau.Value = 10.0;
            SliderCO2Rise.Value = 300.0;
            SliderCO2Breath.Value = 20.0;
            SliderRadarVel.Value = 2.0;
            SliderRadarGain.Value = 1.0;
            SliderRadarClutter.Value = 0.1;
            SliderMicFreq.Value = 200.0;
            SliderMicAmp.Value = 0.5;
            SliderMicHum.Value = 0.05;
            SliderIMUBias.Value = 0.05;
            SliderIMUDrift.Value = 0.01;
            SliderIMUNoise.Value = 0.02;
        }

        private void CalculateAllFFTs()
        {
            var thermal = detectionPipeline.GetSensor<MLX90640Sensor>();
            var co2 = detectionPipeline.GetSensor<CO2Sensor>();
            var radar = detectionPipeline.GetSensor<RadarSensor>();
            var imu = detectionPipeline.GetSensor<IMUSensor>();

            if (sensorSignals.ContainsKey("Thermal")) sensorFFTs["Thermal"] = thermal.CalculateFFT(sensorSignals["Thermal"]);
            if (sensorSignals.ContainsKey("CO2")) sensorFFTs["CO2"] = co2.CalculateFFT(sensorSignals["CO2"]);
            if (sensorSignals.ContainsKey("Radar")) sensorFFTs["Radar"] = radar.CalculateFFT(sensorSignals["Radar"]);
            if (sensorSignals.ContainsKey("Microphone")) sensorFFTs["Microphone"] = CalculateSimpleFFT(sensorSignals["Microphone"]);
            if (sensorSignals.ContainsKey("IMU")) sensorFFTs["IMU"] = imu.CalculateFFT(sensorSignals["IMU"]);
        }

        private Complex[] CalculateSimpleFFT(double[] signal)
        {
            if (signal == null) return Array.Empty<Complex>();
            int n = signal.Length;
            int fftSize = 1;
            while (fftSize < n) fftSize *= 2;

            Complex[] fftData = new Complex[fftSize];
            for (int i = 0; i < n; i++) fftData[i] = new Complex(signal[i], 0);
            for (int i = n; i < fftSize; i++) fftData[i] = Complex.Zero;

            return fftData;
        }

        private void UpdateTimePlot(WpfPlot plot, double[] signal, double fs, string hexColor)
        {
            plot.Plot.Clear();
            if (signal == null || signal.Length == 0) return;

            var sigPlot = plot.Plot.Add.Signal(signal);
            sigPlot.Data.Period = 1.0 / fs;
            sigPlot.Color = ScottPlot.Color.FromHex(hexColor);
            sigPlot.LineWidth = 2f;

            plot.Plot.Axes.AutoScale();
            plot.Refresh();
        }

        private void UpdateFFTPlotDiscrete(WpfPlot plot, Complex[] fft, double fs, string hexColor)
        {
            plot.Plot.Clear();
            if (fft == null || fft.Length == 0) return;

            int N = fft.Length / 2;
            double[] freqs = new double[N];
            double[] mags = new double[N];
            double freqRes = fs / fft.Length;

            for (int i = 0; i < N; i++)
            {
                freqs[i] = i * freqRes;
                double mag = fft[i].Magnitude;
                mags[i] = 20 * Math.Log10(mag + 1e-12);
            }

            double threshold = mags.Max() - 40;

            var color = ScottPlot.Color.FromHex(hexColor);
            for (int i = 0; i < N; i += Math.Max(1, N / 200))
            {
                if (mags[i] > threshold)
                {
                    var line = plot.Plot.Add.Line(freqs[i], -40, freqs[i], mags[i]);
                    line.Color = color;
                    line.LineWidth = 2f;

                    var marker = plot.Plot.Add.Marker(freqs[i], mags[i]);
                    marker.Color = color;
                    marker.Size = 6;
                    marker.Shape = MarkerShape.FilledCircle;
                }
            }

            plot.Plot.Axes.AutoScale();
            plot.Plot.Axes.SetLimitsY(-40, 120);
            plot.Plot.Grid.MajorLineWidth = 1;
            plot.Refresh();
        }

        private void UpdateLaplacePlot()
        {
            PlotLaplace.Plot.Clear();

            // Setup axes
            var vline = PlotLaplace.Plot.Add.VerticalLine(0);
            vline.Color = ScottPlot.Colors.Black;
            vline.LineWidth = 2;

            var hline = PlotLaplace.Plot.Add.HorizontalLine(0);
            hline.Color = ScottPlot.Colors.Black;
            hline.LineWidth = 2;

            // Stability boundary (imaginary axis) - GREEN dashed line
            var stabLine = PlotLaplace.Plot.Add.VerticalLine(0);
            stabLine.Color = ScottPlot.Colors.Green;
            stabLine.LineWidth = 3;
            stabLine.LinePattern = LinePattern.Dashed;

            var allPoles = new List<Complex>();

            try
            {
                var thermal = detectionPipeline.GetSensor<MLX90640Sensor>();
                var co2 = detectionPipeline.GetSensor<CO2Sensor>();
                var radar = detectionPipeline.GetSensor<RadarSensor>();
                var imu = detectionPipeline.GetSensor<IMUSensor>();

                if (thermal != null)
                {
                    var poles = thermal.GetLaplacePolesZeros();
                    allPoles.AddRange(poles);
                }

                if (co2 != null)
                {
                    var poles = co2.GetLaplacePolesZeros();
                    allPoles.AddRange(poles);
                }

                if (radar != null)
                {
                    var poles = radar.GetLaplacePolesZeros();
                    allPoles.AddRange(poles);
                }

                if (imu != null)
                {
                    var poles = imu.GetLaplacePolesZeros();
                    allPoles.AddRange(poles);
                }
            }
            catch (Exception ex)
            {
                LogStatus($"⚠️ Laplace plot error: {ex.Message}");
            }

            // SOLUTION: Use single-point scatter for each pole (NO LINES!)
            if (allPoles.Count > 0)
            {
                LogStatus($"📊 Laplace: Plotting {allPoles.Count} poles");

                foreach (var pole in allPoles)
                {
                    // Create array with ONLY ONE POINT
                    double[] xs = new double[] { pole.Real };
                    double[] ys = new double[] { pole.Imaginary };

                    var scatter = PlotLaplace.Plot.Add.Scatter(xs, ys);
                    scatter.MarkerShape = MarkerShape.Cross;
                    scatter.MarkerSize = 20;
                    scatter.Color = ScottPlot.Colors.Red;
                    scatter.LineWidth = 0;  // NO LINE

                    LogStatus($"  Pole: ({pole.Real:F1}, {pole.Imaginary:F1})");
                }
            }

            // Set proper axis limits
            PlotLaplace.Plot.Axes.SetLimits(-600, 50, -300, 300);

            // Labels
            PlotLaplace.Plot.Axes.Bottom.Label.Text = "Real (σ)";
            PlotLaplace.Plot.Axes.Bottom.Label.FontSize = 12;
            PlotLaplace.Plot.Axes.Left.Label.Text = "Imaginary (jω)";
            PlotLaplace.Plot.Axes.Left.Label.FontSize = 12;
            PlotLaplace.Plot.Title("Laplace S-Plane: X = Poles (Stable if Real < 0)");

            PlotLaplace.Refresh();
        }

        private void UpdateZDomainPlot()
        {
            PlotZDomain.Plot.Clear();

            // Draw unit circle (stability boundary) - GREEN
            var circle = PlotZDomain.Plot.Add.Circle(0, 0, 1);
            circle.LineColor = ScottPlot.Colors.Green;
            circle.LinePattern = LinePattern.Dashed;
            circle.LineWidth = 3;
            circle.FillColor = ScottPlot.Colors.Transparent;

            // Axes
            var vline = PlotZDomain.Plot.Add.VerticalLine(0);
            vline.Color = ScottPlot.Colors.Black;
            vline.LineWidth = 2;

            var hline = PlotZDomain.Plot.Add.HorizontalLine(0);
            hline.Color = ScottPlot.Colors.Black;
            hline.LineWidth = 2;

            var allPoles = new List<Complex>();
            var allZeros = new List<Complex>();

            try
            {
                var thermal = detectionPipeline.GetSensor<MLX90640Sensor>();
                var co2 = detectionPipeline.GetSensor<CO2Sensor>();
                var radar = detectionPipeline.GetSensor<RadarSensor>();
                var imu = detectionPipeline.GetSensor<IMUSensor>();

                if (thermal != null)
                {
                    var pz = thermal.GetZDomainPolesZeros(thermal.SamplingRate);
                    if (pz.Length > 0) allPoles.Add(pz[0]); // pole
                    if (pz.Length > 1 && pz[1].Magnitude > 0.01) allZeros.Add(pz[1]); // zero
                }

                if (co2 != null)
                {
                    var pz = co2.GetZDomainPolesZeros(co2.SamplingRate);
                    if (pz.Length > 0) allPoles.Add(pz[0]);
                    if (pz.Length > 1 && pz[1].Magnitude > 0.01) allZeros.Add(pz[1]);
                }

                if (radar != null)
                {
                    var pz = radar.GetZDomainPolesZeros(radar.SamplingRate);
                    allPoles.AddRange(pz);
                }

                if (imu != null)
                {
                    var pz = imu.GetZDomainPolesZeros(imu.SamplingRate);
                    allPoles.AddRange(pz);
                }
            }
            catch (Exception ex)
            {
                LogStatus($"⚠️ Z-Domain plot error: {ex.Message}");
            }

            // Plot POLES using single-point scatter (X markers in BLUE)
            if (allPoles.Count > 0)
            {
                LogStatus($"📊 Z-Domain: Plotting {allPoles.Count} poles");

                foreach (var pole in allPoles)
                {
                    double[] xs = new double[] { pole.Real };
                    double[] ys = new double[] { pole.Imaginary };

                    var scatter = PlotZDomain.Plot.Add.Scatter(xs, ys);
                    scatter.MarkerShape = MarkerShape.Cross;
                    scatter.MarkerSize = 18;
                    scatter.Color = ScottPlot.Colors.Blue;
                    scatter.LineWidth = 0;

                    double mag = pole.Magnitude;
                    string stable = mag < 1.0 ? "STABLE" : "UNSTABLE";
                    LogStatus($"  Pole: ({pole.Real:F3}, {pole.Imaginary:F3}) |z|={mag:F3} [{stable}]");
                }
            }

            // Plot ZEROS (O markers in RED)
            if (allZeros.Count > 0)
            {
                LogStatus($"📊 Z-Domain: Plotting {allZeros.Count} zeros");

                foreach (var zero in allZeros)
                {
                    double[] xs = new double[] { zero.Real };
                    double[] ys = new double[] { zero.Imaginary };

                    var scatter = PlotZDomain.Plot.Add.Scatter(xs, ys);
                    scatter.MarkerShape = MarkerShape.OpenCircle;
                    scatter.MarkerSize = 18;
                    scatter.Color = ScottPlot.Colors.Red;
                    scatter.LineWidth = 0;

                    LogStatus($"  Zero: ({zero.Real:F3}, {zero.Imaginary:F3})");
                }
            }

            PlotZDomain.Plot.Axes.SquareUnits();
            PlotZDomain.Plot.Axes.SetLimits(-1.5, 1.5, -1.5, 1.5);

            PlotZDomain.Plot.Axes.Bottom.Label.Text = "Real";
            PlotZDomain.Plot.Axes.Bottom.Label.FontSize = 12;
            PlotZDomain.Plot.Axes.Left.Label.Text = "Imaginary";
            PlotZDomain.Plot.Axes.Left.Label.FontSize = 12;
            PlotZDomain.Plot.Title("Z-Domain: X = Poles (stable inside), O = Zeros");

            PlotZDomain.Refresh();
        }

        private void UpdateDetectionStatus(bool detected, double confidence, Dictionary<string, double> confidences)
        {
            if (detected)
            {
                TxtDetectionStatus.Text = "✅ HUMAN DETECTED";
                TxtDetectionStatus.Foreground = new SolidColorBrush(WpfColors.Green);
                BorderDetectionStatus.Background = new SolidColorBrush(WpfColor.FromRgb(232, 245, 233));
            }
            else
            {
                TxtDetectionStatus.Text = "⏹️ NO DETECTION";
                TxtDetectionStatus.Foreground = new SolidColorBrush(WpfColors.Gray);
                BorderDetectionStatus.Background = new SolidColorBrush(WpfColor.FromRgb(250, 250, 250));
            }

            TxtConfidence.Text = $"Overall Confidence: {confidence:F1}%";

            TxtConfThermal.Text = $"🌡️ Thermal: {confidences.GetValueOrDefault("Thermal", 0.0):F1}% " + (confidences.GetValueOrDefault("Thermal", 0.0) >= 50 ? "✓" : "✗");
            TxtConfCO2.Text = $"💨 CO2: {confidences.GetValueOrDefault("CO2", 0.0):F1}% " + (confidences.GetValueOrDefault("CO2", 0.0) >= 50 ? "✓" : "✗");
            TxtConfRadar.Text = $"📡 Radar: {confidences.GetValueOrDefault("Radar", 0.0):F1}% " + (confidences.GetValueOrDefault("Radar", 0.0) >= 50 ? "✓" : "✗");
            TxtConfMic.Text = $"🎤 Mic: {confidences.GetValueOrDefault("Microphone", 0.0):F1}% " + (confidences.GetValueOrDefault("Microphone", 0.0) >= 50 ? "✓" : "✗");
            TxtConfIMU.Text = $"📊 IMU: {confidences.GetValueOrDefault("IMU", 0.0):F1}% " + (confidences.GetValueOrDefault("IMU", 0.0) >= 50 ? "✓" : "✗");

            LogStatus($"Detection: {(detected ? "POSITIVE" : "NEGATIVE")} | Confidence: {confidence:F1}%");
        }

        private void UpdateSignalMetrics()
        {
            if (!sensorSignals.ContainsKey("Thermal") || sensorSignals["Thermal"].Length == 0) return;

            var thermal = detectionPipeline.GetSensor<MLX90640Sensor>();
            var signal = sensorSignals["Thermal"];

            TxtMetricRMS.Text = $"RMS: {thermal.CalculateRMS(signal):F3}";
            TxtMetricPeakToPeak.Text = $"Peak-to-Peak: {thermal.CalculatePeakToPeak(signal):F3}";
            TxtMetricMean.Text = $"Mean: {thermal.CalculateMean(signal):F3}";
            TxtMetricStdDev.Text = $"Std Dev: {thermal.CalculateStdDev(signal):F3}";
        }

        private void LogStatus(string message)
        {
            Dispatcher.Invoke(() =>
            {
                string timestamp = DateTime.Now.ToString("HH:mm:ss");
                TxtStatusLog.Text += $"\n[{timestamp}] {message}";
                StatusScroller.ScrollToEnd();
            });
        }
    }
}