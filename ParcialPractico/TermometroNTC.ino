#include <WiFi.h>
#include <ThingSpeak.h>
#include <Adafruit_ADS1X15.h>
#include "AdafruitIO_WiFi.h"

// ----------------- CONFIG (NO subir a repos publicos) -----------------
#define IO_USERNAME  "Nicolas804"
#define IO_KEY       "aio_GxCR21RprJ6SKYBAHcde3SR90fAj"
#define WIFI_SSID    "FAMILIA LESMES"
#define WIFI_PASS    "hugoyandrea2021"

constexpr unsigned long CHANNEL_ID = 3167257UL;
constexpr const char* WRITE_API_KEY = "CLPWTX8QWH65HFIL";

// ----------------- HARDWARE / LIBS -----------------
Adafruit_ADS1115 ads;
WiFiClient client;

// Adafruit IO
AdafruitIO_WiFi io(IO_USERNAME, IO_KEY, WIFI_SSID, WIFI_PASS);
AdafruitIO_Feed *feed_temp = io.feed("TemperaturaMedicion");

// ADS1115 LSB en mV (PGA por defecto  ±6.144V -> 0.1875 mV/bit)
constexpr float LSB_MV = 0.1875f;         // mV por bit
constexpr float MV_TO_V = 0.001f;        // para convertir mV a V

// ------------- FUNCIONES -------------
float calcularTemperatura(float volts) {
  // volts en V
  if (volts >= 4.93f) {
    return (-1.0f / 0.035354f) * logf((5.007271f - volts) / 1.691408f);
  } else {
    return (logf(1.0f / (1.0f - ((volts - 0.293166f) / 4.665771f))) / 0.05133f) - 7.9f;
  }
}

void conectarWiFiConTimeout(unsigned long timeoutMs = 20000UL) {
  WiFi.mode(WIFI_STA);
  WiFi.disconnect(true);
  delay(100);

  Serial.printf("Conectando a '%s' ...\n", WIFI_SSID);
  WiFi.begin(WIFI_SSID, WIFI_PASS);
  unsigned long start = millis();
  while (WiFi.status() != WL_CONNECTED && (millis() - start) < timeoutMs) {
    Serial.print('.');
    delay(300);
  }
  Serial.println();
  if (WiFi.status() == WL_CONNECTED) {
    Serial.print(F("WiFi conectado, IP: "));
    Serial.println(WiFi.localIP());
  } else {
    Serial.println(F("No conectado a WiFi (timeout)."));
  }
}

void setup() {
  Serial.begin(115200);
  delay(300);

  Serial.println(F("Iniciando..."));
  conectarWiFiConTimeout(20000UL);

  ThingSpeak.begin(client);

  // Conectar Adafruit IO (no bloqueante excesivo)
  Serial.println(F("Conectando a Adafruit IO..."));
  io.connect();
  unsigned long start = millis();
  const unsigned long aioTimeout = 15000UL;
  while (io.status() < AIO_CONNECTED && (millis() - start) < aioTimeout) {
    Serial.print('.');
    io.run(); // mantener el proceso de conexión
    delay(300);
  }
  Serial.println();
  if (io.status() >= AIO_CONNECTED) {
    Serial.println(F("Adafruit IO conectado."));
  } else {
    Serial.println(F("No se pudo conectar a Adafruit IO en el tiempo dado."));
  }

  // I2C y ADS1115
  Wire.begin(21, 22);
  if (!ads.begin()) {
    Serial.println(F("Error inicializando ADS1115."));
    while (true) { delay(10); }
  }
  Serial.println(F("ADS1115 listo."));
}

void loop() {
  // Leer diferencial 0-1
  int16_t raw = ads.readADC_Differential_0_1();

  // Convertir a voltios: raw * LSB_MV (mV) -> V
  float volts = (raw * LSB_MV) * MV_TO_V;

  // Calcular temperatura (solo una vez)
  float tempC = calcularTemperatura(volts);

  // Mostrar info (evitar Strings)
  Serial.printf("Raw: %d  Voltios: %.4f V  Temp: %.3f C\n", raw, volts, tempC);

  // Enviar a ThingSpeak si hay WiFi
  if (WiFi.status() == WL_CONNECTED) {
    ThingSpeak.setField(1, tempC);
    int tsResult = ThingSpeak.writeFields(CHANNEL_ID, WRITE_API_KEY);
    if (tsResult == 200) {
      Serial.println(F("ThingSpeak: OK"));
    } else {
      Serial.printf("ThingSpeak error: %d\n", tsResult);
    }
  } else {
    Serial.println(F("WiFi no conectado -> no se envía a ThingSpeak."));
  }

  // Enviar a Adafruit IO si está conectado
  // `io.run()` mantiene la conexión y procesa callbacks; debe llamarse frecuentemente.
  io.run();
  if (io.status() >= AIO_CONNECTED) {
    feed_temp->save(tempC);
    Serial.println(F("Adafruit IO: dato enviado (feed save)."));
  } else {
    Serial.println(F("Adafruit IO no conectado -> no se envía al feed."));
  }

  // Pequeña pausa
  delay(2500);
}
