// Arduino code for Ardunio Nano 33 BLE Sense 
// This code sets up the Arduino to work with a SSD1306 Display
// and a stepper motor controller. The arduino will accept commands
// over serial and perform the following actions
//
// I -> Start displaying IMU data
// i -> Stop displaying IMU data
// C -> Start displaying color data
// c -> Stop displaying color data
// G -> Spin stepper motor according to current setting of numSteps and stepTime, where
//      numSteps are the number of steps, and stepTime is the time in micro seconds per step
// R -> Set for turn CW (Right), needs to be passed 2 two byte integer which specify
//      the numSteps and stepTime. e.g. Rxxyy (where xx and yy are little endian integers)
// r -> Set for turn CW (Right), same as R, but uses ascii digits encoded as 
//      rnnnn,nnnn; where the first number is numSteps and the second number is stepTime
//      e.g. r800,500; is right turn for 800 steps with 500 microseconds per step
// L,l -> same as R,r but for CCW (Left) turns

////////////////////
// Include section
////////////////////
#include <Adafruit_SSD1306.h>
#include <Arduino_BMI270_BMM150.h>
#include <Arduino_APDS9960.h>

/////////////////////
// global variables
/////////////////////
bool imuEnable = false;
bool colorSensorEnable = false;
bool runStepperMotor = false;
bool displayUpdate = false;

////////////////////////////////////////////////////////////////
// setup
////////////////////////////////////////////////////////////////
void setup() {
  
  // Start serial communication
  setupSerial();

  // setup OLED display (see function below for parameter settings)
  setupDisplay();

  // setup IMU
  setupIMU();

  //setup StepperMotor
  setupStepperMotor();

  // setup light sensor
  setupColorSensor();
}

////////////////////////////////////////////////////////////////
// loop function
////////////////////////////////////////////////////////////////
void loop() {
  // check serial for commands
  updateSerial();

  // display light sensor data
  if (colorSensorEnable) updateColorSensor();

  // display imu if called for
  if (imuEnable) displayIMU();
  
  // run stepper motor if called for
  if (runStepperMotor) spinStepperMotor();

  // flush the display if necessary
  if (displayUpdate) displayFlush();
  displayUpdate = false;
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

///////////////////
// display IMU
///////////////////
void displayIMU() {
  // x, y, z values for sensor info
  float x, y, z;

  if (!displayUpdate) {
    // clear display
    displayClear();
  
    // sample rate display
    displayPrintln(sampleRatesDisplayBuffer);

    // set display update
    displayUpdate = true;
  }

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
}

//////////////////////
// setupStepperMotor
//////////////////////
#define dirPin 2
#define stepPin 3
#define enablePin 4

// some variables that can be set to control the spin
int stepDir = 1;
int numSteps = 800;
int stepTime = 500;

////////////////////////////
// setup the stepper motor
////////////////////////////
void setupStepperMotor() {
  // Sets the enable pin as output and disable the motor
  pinMode(enablePin,OUTPUT);
  digitalWrite(enablePin,HIGH);
  // setup step and dir pins
  pinMode(stepPin,OUTPUT); 
  pinMode(dirPin,OUTPUT);
  // enable the motor
  digitalWrite(enablePin,LOW);
}

////////////////////////////
// spin the stepper motor
////////////////////////////
void spinStepperMotor() {
    // Enables the motor to move in a particular direction
  digitalWrite(dirPin, stepDir == 1 ? HIGH : LOW); 
  
  // keep time for how long this ran
  unsigned long startTime = micros();

  // Run for set number of steps
  for(int iStep = 0; iStep < numSteps; iStep++) {
    digitalWrite(stepPin,HIGH); 
    delayMicroseconds(stepTime/2);    // by changing this time delay between the steps we can change the rotation speed
    digitalWrite(stepPin,LOW); 
    delayMicroseconds(stepTime/2); 
  }
  unsigned long elapsedTime = micros() - startTime;
  
  // clear display and show what we have done
  displayClear();
  displayUpdate = true;
  displayPrintln(dirPin == 1 ? "Right" : "Left");
  displayPrintln("numSteps: " + String(numSteps));
  displayPrintln("stepTime: " + String(stepTime) + " us");
  displayPrintln("Dur: " + String(elapsedTime/1000) + " ms");

  
  // stepper motor done running
  runStepperMotor = false;
}

////////////////////////////////
// Serial Communications
////////////////////////////////
int baudRate = 9600;

////////////////
// setupSerial
////////////////
void setupSerial() {
  // set baud rate
  Serial.begin(baudRate);
  // print that we started
  Serial.println("Init");
}

/////////////////
// udpateSerial
/////////////////
void updateSerial() {
  // check if there are any bytes available
  if (Serial.available() > 0) { 
    
    // clear screen
    displayUpdate = true;
    displayClear();

    char command = Serial.read(); 
    switch(command) {
      case 'R':
      case 'L':
        // set the direction
        stepDir = toupper(command) == 'L' ? 1 : -1;
        // read one int for how many steps
        numSteps = readIntSerial();
        // read another int for the cycle time
        stepTime = readIntSerial();
        // display what we got
        displayPrintln("Received R");
        displayPrintln("numSteps: "+String(numSteps));
        displayPrintln("stepTime: "+String(stepTime));
        break;
    case 'r':
    case 'l':
        // lower case signifies getting the values by ascii
        // so the values you should be given as nnn,nnn;
        // set the direction
        stepDir = toupper(command) == 'L' ? 1 : -1;
        // read one int for how many steps
        numSteps = readIntAsciiSerial(',');
        // read another int for the cycle time
        stepTime = readIntAsciiSerial(';');
        // display what we got
        displayPrintln("Received R");
        displayPrintln("numSteps: "+String(numSteps));
        displayPrintln("stepTime: "+String(stepTime));
        break;
      // GO stepper motor
      case 'G':
        runStepperMotor = true;
        break;
      // IMU 
      case 'I':
        displayPrintln("IMU enabled");
        imuEnable = true;
        break;
      case 'i':
        displayPrintln("IMU disabled");
        imuEnable = false;
        break;
      // Color Sensor
      case 'C':
        displayPrintln("Color Sensor enable");
        colorSensorEnable = true;
        break;
      case 'c':
        displayPrintln("Color Sensor disabled");
        colorSensorEnable = false;
        break;
      default:
        displayPrintln("Unrecogonized");
        break;

    }

  }
}

//////////////////
// readIntSerial
//////////////////
int readIntSerial() {
  // wait for two bytes to become available
  while (Serial.available() < 2) {}
  
  // read the two bytes
  byte byte1 = Serial.read(); // Read the first byte
  byte byte2 = Serial.read(); // Read the second byte

  // turn into an int
  int receivedInt = (byte2 << 8) | byte1;
  return(receivedInt);
}

//////////////////
// readIntAsciiSerial
//////////////////
int readIntAsciiSerial(char delimitter) {
  // read until delimitter
  String receivedString = Serial.readStringUntil(delimitter);
  
  // convert to int
  return(receivedString.toInt());
}

/////////////////////
// setupColorSensor
/////////////////////
void setupColorSensor() {
 if (!APDS.begin()) {
    displayString("Error initializing APDS-9960 sensor.");
  }
}

/////////////////////
// updateColorSensor
/////////////////////
String colorSensorString;
void updateColorSensor() {
  // see whether screen needs to be cleared
  if (!displayUpdate) {
    displayClear();
    displayUpdate = true;
  }
  
  // check for color
  if (APDS.colorAvailable()) {
    int r, g, b, a;
    APDS.readColor(r, g, b, a);
    // display the values
    colorSensorString = String(r)+", "+String(g)+", "+String(b)+", "+String(a);
  }
  // display light sensor string
  displayPrintln(colorSensorString);

}