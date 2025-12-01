using System;
using System.Linq;
using System.Collections.Generic;
using MultiSensorSimulation.Models;

namespace MultiSensorSimulation.Utils
{
    public class SensorFusion
    {
        // Weighted confidence calculation untuk masing-masing sensor
        private double thermalWeight = 0.25;
        private double co2Weight = 0.20;
        private double radarWeight = 0.30;
        private double audioWeight = 0.15;
        private double motionWeight = 0.10;

        // Threshold untuk deteksi
        public double ThermalThreshold { get; set; } = 32.0;  // °C (di atas ambient)
        public double CO2Threshold { get; set; } = 800.0;      // ppm (elevated)
        public double RadarVarianceThreshold { get; set; } = 0.08;
        public double AudioRMSThreshold { get; set; } = 0.12;
        public double MotionPeakThreshold { get; set; } = 0.8;

        /// <summary>
        /// Menghitung confidence dari thermal sensor berdasarkan temperature signature
        /// </summary>
        public double CalculateThermalConfidence(double[] thermalSignal)
        {
            if (thermalSignal == null || thermalSignal.Length == 0)
                return 0.0;

            double avgTemp = thermalSignal.Average();
            double maxTemp = thermalSignal.Max();

            // Confidence meningkat saat mendekati body temperature
            double tempDiff = maxTemp - ThermalThreshold;
            double confidence = Math.Min(100.0, Math.Max(0.0, (tempDiff / 5.0) * 100.0));

            return confidence;
        }

        /// <summary>
        /// Menghitung confidence dari CO2 sensor berdasarkan level concentration
        /// </summary>
        public double CalculateCO2Confidence(double[] co2Signal)
        {
            if (co2Signal == null || co2Signal.Length == 0)
                return 0.0;

            double avgCO2 = co2Signal.Average();
            double maxCO2 = co2Signal.Max();

            // Breathing pattern detection (variability)
            double variance = CalculateVariance(co2Signal);
            double breathingScore = Math.Min(1.0, variance / 100.0);

            // Level confidence
            double levelDiff = maxCO2 - CO2Threshold;
            double levelScore = Math.Min(1.0, Math.Max(0.0, levelDiff / 300.0));

            double confidence = ((levelScore * 0.7) + (breathingScore * 0.3)) * 100.0;
            return Math.Min(100.0, confidence);
        }

        /// <summary>
        /// Menghitung confidence dari radar berdasarkan Doppler signature
        /// </summary>
        public double CalculateRadarConfidence(double[] radarSignal)
        {
            if (radarSignal == null || radarSignal.Length == 0)
                return 0.0;

            // Movement variance detection
            double variance = CalculateVariance(radarSignal);
            double peakToPeak = radarSignal.Max() - radarSignal.Min();

            // Doppler signature strength
            double varianceScore = Math.Min(1.0, variance / RadarVarianceThreshold);
            double amplitudeScore = Math.Min(1.0, peakToPeak / 2.0);

            double confidence = ((varianceScore * 0.6) + (amplitudeScore * 0.4)) * 100.0;
            return Math.Min(100.0, confidence);
        }

        /// <summary>
        /// Menghitung confidence dari microphone berdasarkan audio signature
        /// </summary>
        public double CalculateAudioConfidence(double[] audioSignal)
        {
            if (audioSignal == null || audioSignal.Length == 0)
                return 0.0;

            double rms = CalculateRMS(audioSignal);
            double peakToPeak = audioSignal.Max() - audioSignal.Min();

            // Voice/breathing sound detection
            double rmsScore = Math.Min(1.0, rms / AudioRMSThreshold);
            double dynamicScore = Math.Min(1.0, peakToPeak / 1.0);

            double confidence = ((rmsScore * 0.7) + (dynamicScore * 0.3)) * 100.0;
            return Math.Min(100.0, confidence);
        }

        /// <summary>
        /// Menghitung confidence dari IMU berdasarkan motion signature
        /// </summary>
        public double CalculateMotionConfidence(double[] imuSignal)
        {
            if (imuSignal == null || imuSignal.Length == 0)
                return 0.0;

            // Extract accelerometer data (even indices)
            double[] accelData = imuSignal.Where((x, i) => i % 2 == 0).ToArray();

            double variance = CalculateVariance(accelData);
            double peakToPeak = accelData.Max() - accelData.Min();

            double varianceScore = Math.Min(1.0, variance / 0.5);
            double amplitudeScore = Math.Min(1.0, peakToPeak / MotionPeakThreshold);

            double confidence = ((varianceScore * 0.5) + (amplitudeScore * 0.5)) * 100.0;
            return Math.Min(100.0, confidence);
        }

        /// <summary>
        /// Fusion akhir dengan weighted sum semua sensor confidences
        /// </summary>
        public double CalculateFinalConfidence(
            double c_thermal,
            double c_co2,
            double c_radar,
            double c_audio,
            double c_motion)
        {
            double finalConfidence =
                (thermalWeight * c_thermal) +
                (co2Weight * c_co2) +
                (radarWeight * c_radar) +
                (audioWeight * c_audio) +
                (motionWeight * c_motion);

            return Math.Min(100.0, Math.Max(0.0, finalConfidence));
        }

        /// <summary>
        /// Deteksi victim berdasarkan final confidence dan jumlah sensor positif
        /// </summary>
        public bool IsVictimDetected(
            double finalConfidence,
            int positiveSensors,
            double confidenceThreshold = 70.0,
            int minPositiveSensors = 3)
        {
            return finalConfidence >= confidenceThreshold &&
                   positiveSensors >= minPositiveSensors;
        }

        /// <summary>
        /// Hitung berapa banyak sensor yang memberikan positive detection
        /// </summary>
        public int CountPositiveSensors(
            double c_thermal,
            double c_co2,
            double c_radar,
            double c_audio,
            double c_motion,
            double threshold = 50.0)
        {
            int count = 0;
            if (c_thermal >= threshold) count++;
            if (c_co2 >= threshold) count++;
            if (c_radar >= threshold) count++;
            if (c_audio >= threshold) count++;
            if (c_motion >= threshold) count++;
            return count;
        }

        // Helper methods
        private double CalculateVariance(double[] data)
        {
            double mean = data.Average();
            return data.Select(x => Math.Pow(x - mean, 2)).Average();
        }

        private double CalculateRMS(double[] data)
        {
            return Math.Sqrt(data.Select(x => x * x).Average());
        }
    }
}