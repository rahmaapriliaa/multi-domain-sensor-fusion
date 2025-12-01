using System;
using System.Collections.Generic;
using System.Linq;
using MultiSensorSimulation.Models;
using MultiSensorSimulation.Utils;

namespace MultiSensorSimulation.Utils
{
    public class DetectionPipeline
    {
        // ========================================
        // SENSOR INSTANCES
        // ========================================
        // FIX: Tambahkan = null! untuk mengatasi error CS8618
        private MLX90640Sensor thermalSensor = null!;
        private CO2Sensor co2Sensor = null!;
        private RadarSensor radarSensor = null!;
        private MicSensor micSensor = null!;
        private IMUSensor imuSensor = null!;

        // Sensor Fusion Engine
        private SensorFusion fusionEngine = null!;

        // ========================================
        // CONFIGURABLE PARAMETERS
        // ========================================
        public double ConfidenceThreshold { get; set; } = 70.0;
        public int MinimumPositiveSensors { get; set; } = 3;

        public bool EnableThermal { get; set; } = true;
        public bool EnableCO2 { get; set; } = true;
        public bool EnableRadar { get; set; } = true;
        public bool EnableMicrophone { get; set; } = true;
        public bool EnableIMU { get; set; } = true;

        // ========================================
        // CONSTRUCTOR
        // ========================================
        public DetectionPipeline()
        {
            fusionEngine = new SensorFusion();
            InitializeSensors();
        }

        private void InitializeSensors()
        {
            thermalSensor = new MLX90640Sensor
            {
                Name = "MLX90640 Thermal",
                Unit = "°C",
                Emissivity = 0.95,
                TimeConstant = 2.0,
                AmbientTemp = 25.0,
                BodyTemp = 37.0,
                SamplingRate = 16.0
            };

            co2Sensor = new CO2Sensor
            {
                Name = "CO2 Sensor",
                Unit = "ppm",
                Baseline = 450.0,
                HumanRise = 300.0,
                TimeConstantResp = 10.0,
                SamplingRate = 10.0
            };

            radarSensor = new RadarSensor
            {
                Name = "24GHz Radar",
                Unit = "V",
                CarrierFrequency = 24.0e9,
                TargetVelocity = 2.0,
                SamplingRate = 1000.0
            };

            // Pastikan MicSensor memiliki properti ini (lihat Solusi 3)
            micSensor = new MicSensor
            {
                Name = "MEMS Microphone",
                Unit = "V",
                SignalAmplitude = 0.5,
                SignalFrequency = 200.0,
                HumFrequency = 50.0
            };

            imuSensor = new IMUSensor
            {
                Name = "6-Axis IMU",
                Unit = "m/s²",
                AccelSensitivity = 1.0,
                GyroSensitivity = 1.0,
                SamplingRate = 100.0
            };
        }

        public Dictionary<string, double[]> RunAllSensors(double gain, double tau, double duration)
        {
            var results = new Dictionary<string, double[]>();

            try
            {
                if (EnableThermal) results["Thermal"] = thermalSensor.GenerateSignal(gain, tau, duration);
                if (EnableCO2) results["CO2"] = co2Sensor.GenerateSignal(gain, tau, duration);
                if (EnableRadar) results["Radar"] = radarSensor.GenerateSignal(gain, tau, duration);
                if (EnableMicrophone) results["Microphone"] = micSensor.GenerateSignal(gain, tau, duration);
                if (EnableIMU) results["IMU"] = imuSensor.GenerateSignal(gain, tau, duration);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error in sensor acquisition: {ex.Message}");
                throw;
            }

            return results;
        }

        public T GetSensor<T>() where T : SensorBase
        {
            if (typeof(T) == typeof(MLX90640Sensor)) return thermalSensor as T;
            if (typeof(T) == typeof(CO2Sensor)) return co2Sensor as T;
            if (typeof(T) == typeof(RadarSensor)) return radarSensor as T;
            if (typeof(T) == typeof(MicSensor)) return micSensor as T;
            if (typeof(T) == typeof(IMUSensor)) return imuSensor as T;
            return null;
        }

        public (bool detected, double confidence, Dictionary<string, double> confidences) DetectHuman(Dictionary<string, double[]> sensorData)
        {
            var confidences = CalculateIndividualConfidences(sensorData);

            double c_thermal = confidences.GetValueOrDefault("Thermal", 0.0);
            double c_co2 = confidences.GetValueOrDefault("CO2", 0.0);
            double c_radar = confidences.GetValueOrDefault("Radar", 0.0);
            double c_audio = confidences.GetValueOrDefault("Microphone", 0.0);
            double c_motion = confidences.GetValueOrDefault("IMU", 0.0);

            double finalConfidence = fusionEngine.CalculateFinalConfidence(c_thermal, c_co2, c_radar, c_audio, c_motion);
            int positiveSensors = fusionEngine.CountPositiveSensors(c_thermal, c_co2, c_radar, c_audio, c_motion, threshold: 50.0);

            bool detected = fusionEngine.IsVictimDetected(finalConfidence, positiveSensors, ConfidenceThreshold, MinimumPositiveSensors);

            return (detected, finalConfidence, confidences);
        }

        public Dictionary<string, double> CalculateIndividualConfidences(Dictionary<string, double[]> sensorData)
        {
            var confidences = new Dictionary<string, double>();

            if (sensorData.ContainsKey("Thermal")) confidences["Thermal"] = fusionEngine.CalculateThermalConfidence(sensorData["Thermal"]);
            if (sensorData.ContainsKey("CO2")) confidences["CO2"] = fusionEngine.CalculateCO2Confidence(sensorData["CO2"]);
            if (sensorData.ContainsKey("Radar")) confidences["Radar"] = fusionEngine.CalculateRadarConfidence(sensorData["Radar"]);
            if (sensorData.ContainsKey("Microphone")) confidences["Microphone"] = fusionEngine.CalculateAudioConfidence(sensorData["Microphone"]);
            if (sensorData.ContainsKey("IMU")) confidences["IMU"] = fusionEngine.CalculateMotionConfidence(sensorData["IMU"]);

            return confidences;
        }
    }
}