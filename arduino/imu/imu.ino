#include <Adafruit_SSD1306.h>
#include "Arduino_BMI270_BMM150.h"

//#include <splash.h>
//#include <Nano33BLE_System.h>

////////////////////////////////////////////////////////////////
// setup
////////////////////////////////////////////////////////////////
void setup() {
  
  // Start serial communication
  Serial.begin(9600);
  Serial.println("Initializing program");

  // setup OLED display (see function below for parameter settings)
  setupDisplay();

  // setup IMU
  //setupIMU();

  //setup StepperMotor
  setupStepperMotor();
}

////////////////////////////////////////////////////////////////
// loop function
////////////////////////////////////////////////////////////////
void loop() {
  //displayIMU();
  spinStepperMotor();
}

//////////////////////////////////////////////////////
// OLED display
//////////////////////////////////////////////////////
// parameters
#define SCREEN_WIDTH 128     // OLED display width, in pixels
#define SCREEN_HEIGHT 32     // OLED display height, in pixels
#define OLED_RESET -1        // Reset pin # (or -1 if sharing Arduino reset pin)
#define TOP_SCREEN_ADDRESS 0x3D  ///< See datasheet for Address; 0x3D for 128x64, 0x3C for 128x32
#define BOTTOM_SCREEN_ADDRESS 0x3C  ///< See datasheet for Address; 0x3D for 128x64, 0x3C for 128x32
#define TEXT_HEIGHT 8

int nDisplays = 0;
Adafruit_SSD1306 displays[] = {
  Adafruit_SSD1306(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, OLED_RESET),
  Adafruit_SSD1306(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, OLED_RESET)
};

///////////////////
// setup function
///////////////////
void setupDisplay() {
  // SSD1306_SWITCHCAPVCC = generate display voltage from 3.3V internally
  if (!displays[nDisplays].begin(SSD1306_SWITCHCAPVCC, TOP_SCREEN_ADDRESS)) {
    Serial.println(F("(setupDisplay) Could not initialized top display: SSD1306 allocation failed"));
  }
  else {
    displays[nDisplays].setRotation(2);
    nDisplays++;
  }

  if (!displays[nDisplays].begin(SSD1306_SWITCHCAPVCC, BOTTOM_SCREEN_ADDRESS)) {
    Serial.println(F("(setupDisplay) Could not initialized bottom display: SSD1306 allocation failed"));
  }
  else {
    nDisplays++;
  }

  // Display Hello World Message
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    displays[iDisplay].clearDisplay();
    displays[iDisplay].setTextSize(1);               // Normal 1:1 pixel scale
    displays[iDisplay].setTextColor(SSD1306_WHITE);  // Draw white text
    displays[iDisplay].setCursor(0, 0);              // Start at top-left corner
    displays[iDisplay].println("Hello World");
    displays[iDisplay].display();
  }
}

/////////////////
// flush display
/////////////////
void displayFlush() {
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    // flush
    displays[iDisplay].display();
  }
}

/////////////////
// cleardisplay
/////////////////
void displayClear() {
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    displays[iDisplay].clearDisplay();
    displays[iDisplay].setCursor(0,0);
  }
}
/////////////////
// display a string
/////////////////
void displayString(String str) {
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    displays[iDisplay].clearDisplay();
    displays[iDisplay].setCursor(0, 0);              // Start at top-left corner
    displays[iDisplay].println(str);
    displays[iDisplay].display();
  }
}

/////////////////
// displayPrint
/////////////////
void displayPrint(String str) {
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    displays[iDisplay].print(str);
  }
}

///////////////////
// displayPrintln
///////////////////
void displayPrintln(String str) {
  for (int iDisplay = 0; iDisplay<nDisplays; iDisplay++) {
    displays[iDisplay].println(str);
  }
}

///////////////////
// setup IMU
///////////////////
// buffers for display
char sampleRatesDisplayBuffer[32];
char magDisplayBuffer[32];
char accDisplayBuffer[32];
char gyroDisplayBuffer[32];

void setupIMU() {
  if (!IMU.begin()) {
    displayPrintln("Failed to initialize IMU!");
    displayFlush();
    while (1);
  }
  // get sample rates
  sprintf(sampleRatesDisplayBuffer, "M:%04.1f G:%04.1f A:%05.2f", IMU.magneticFieldSampleRate(), IMU.gyroscopeSampleRate(), IMU.accelerationSampleRate());
}

void displayIMU() {
  // x, y, z values for sensor info
  float x, y, z;

  // clear display
  displayClear();
  
  // sample rate display
  displayPrintln(sampleRatesDisplayBuffer);

  // print magnetic field
  if (IMU.magneticFieldAvailable()) {
    // read magnetic field
    IMU.readMagneticField(x, y, z);
    // create display string
    sprintf(magDisplayBuffer, "M:%+05.1f %+05.1f %+05.1f", x, y, z);
  }
  displayPrintln(magDisplayBuffer);

  // print gyroscope
  if (IMU.gyroscopeAvailable()) {
    // read magnetic field
    IMU.readGyroscope(x, y, z);
    // display
    sprintf(gyroDisplayBuffer, "G:%+05.1f %+05.1f %+05.1f", x, y, z);
  }
  displayPrintln(gyroDisplayBuffer);

  // print accelerometer
  if (IMU.accelerationAvailable()) {
    // read magnetic field
    IMU.readAcceleration(x, y, z);
    // display
    sprintf(accDisplayBuffer, "A:%+05.2f %+05.2f %+05.2f", x, y, z);
  }
  displayPrintln(accDisplayBuffer);
  displayFlush();
}

//////////////////////
// setupStepperMotor
//////////////////////
#define stepPin 3
#define dirPin 2
#define enPin 4
void setupStepperMotor() {
  // Sets the two pins as Outputs
  pinMode(enPin,OUTPUT);
  digitalWrite(enPin,HIGH);
  pinMode(stepPin,OUTPUT); 
  pinMode(dirPin,OUTPUT);
  digitalWrite(enPin,LOW);
}

void spinStepperMotor() {
  digitalWrite(dirPin,HIGH); // Enables the motor to move in a particular direction
  displayString("Spin");
  // Makes 200 pulses for making one full cycle rotation
  for(int x = 0; x < 800; x++) {
    digitalWrite(stepPin,HIGH); 
    delayMicroseconds(700);    // by changing this time delay between the steps we can change the rotation speed
    digitalWrite(stepPin,LOW); 
    delayMicroseconds(700); 
  }
  delay(1000); // One second delay
  displayString("Reverse Spin");
  
  digitalWrite(dirPin,LOW); //Changes the rotations direction
  // Makes 400 pulses for making two full cycle rotation
  for(int x = 0; x < 1600; x++) {
    digitalWrite(stepPin,HIGH);
    delayMicroseconds(500);
    digitalWrite(stepPin,LOW);
    delayMicroseconds(500);
  }
  delay(1000);
}