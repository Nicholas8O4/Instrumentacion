#include <Wire.h>
#include <Adafruit_ADS1X15.h>

Adafruit_ADS1115 ads;

// ---------------- FUNCIONES ----------------

// Convierte un voltaje a lux y devuelve un mensaje formateado
String calcularLuxFotorresistencia(float volts) {
  float lux =  (3717.988649/(5.481759-volts))-281.186480;  

  // Crear el mensaje
  String mensaje = "Luminosidad (LUX): ";
  mensaje += String(lux, 2);   // 2 decimales es suficiente para lux
  mensaje += " lx";

  return mensaje;
}


String calcularLuxFotodiodo(float volts) {
  float lux =  (49968.732877/(15.608029-volts))-3217.375803;  

  // Crear el mensaje
  String mensaje = "Luminosidad (LUX): ";
  mensaje += String(lux, 2);   // 2 decimales es suficiente para lux
  mensaje += " lx";

  return mensaje;
}


String calcularPesoGalgaResistiva(float volts) {
  float peso = (732.528382/(4.787676-volts))+39.040351;  

  // Crear el mensaje
  String mensaje = "Peso (g): ";
  mensaje += String(peso, 2);   // 2 decimales es suficiente para lux
  mensaje += " g";

  return mensaje;
}

// ---------------- SETUP ----------------

void setup() {
  Serial.begin(115200);
  delay(50);

  Wire.begin(21, 22); // Pines I2C para ESP32

  Serial.println("Inicializando ADS1115...");

  if (!ads.begin()) {
    Serial.println("Error: No se pudo inicializar el ADS1115.");
    while (1) { delay(10); }
  }

  Serial.println("ADS1115 listo para medir.");
  Serial.println("-----------------------------------------------------------");
}

// ---------------- LOOP ----------------

void loop() {
  // Lecturas de los canales A0, A1 y A2
  int16_t adc0 = ads.readADC_SingleEnded(0);
  int16_t adc1 = ads.readADC_SingleEnded(1);
  int16_t adc2 = ads.readADC_SingleEnded(2);

  // Conversión a voltios
  float volts0 = ads.computeVolts(adc0);
  float volts1 = ads.computeVolts(adc1);
  float volts2 = ads.computeVolts(adc2);

  float volts00 = (adc0 * 6.144) / 32767.0;

  // Impresión de valores crudos y convertidos
  Serial.println("Lecturas actuales:");
  Serial.print("AIN0: "); Serial.print(adc0); Serial.print("\t"); Serial.print(volts00, 4); Serial.println(" V");
  Serial.print("AIN1: "); Serial.print(adc1); Serial.print("\t"); Serial.print(volts1, 4); Serial.println(" V");
  Serial.print("AIN2: "); Serial.print(adc2); Serial.print("\t"); Serial.print(volts2, 4); Serial.println(" V");

  Serial.println("-----------------------------------------------------------");
  // Evaluar cada canal con la función calcularLux()
  Serial.println(calcularLuxFotorresistencia(volts0));
  Serial.println(calcularPesoGalgaResistiva(volts1));
  Serial.println(calcularLuxFotodiodo(volts2));

  Serial.println("-----------------------------------------------------------");
  delay(1000);
}
