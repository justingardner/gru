///////////////
// constants //
///////////////

// This pin goes to the MOSFET to control the LED
int ledControlPin = 7;

// This pin goes to the test button
// and has a pull-down resistor so
// that is off by default
int buttonPin = 3;

// this pin is connected to the BNC port
// and is for detecting triggers. Note
// that it uses an interrupt, so can only
// be pin or 2 or 3 on an Uno board
int triggerPin = 2;

// gTriggerMillis is the time at
// which a trigger was encountered
// and is set by the triggerInterrupt function
volatile unsigned long gTriggerMillis = 0;
// gTrigger is whether a trigger has
// been encountered
volatile int gTrigger = false;

// stimulus control
int iCycle;
int nCycles = 5;
int ledOnTime = 1;
int ledOffTime = 100;

////////////////////
// setup function //
////////////////////
void setup() {
  // initialize ledControlPin to be an output pin
  // and set to low, so that LED starts out off
  pinMode(ledControlPin, OUTPUT);
  digitalWrite(ledControlPin, LOW);

  // attach the interrupt to detect digital pulses
  // (the rising phase) on the triggerPin and buttonPin
  // then call the interrupt function, triggerInterrupt
  attachInterrupt(digitalPinToInterrupt(triggerPin), triggerInterrupt, RISING);
  attachInterrupt(digitalPinToInterrupt(buttonPin), triggerInterrupt, RISING);

  // set the builtin LED to output
  // and turn off the LED
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(LED_BUILTIN, LOW);
}

///////////////////
// loop function //
///////////////////
void loop() {

  if (gTrigger) {
    //Turn built in LED on
    digitalWrite(LED_BUILTIN, HIGH);

    // Cycle for LED display
    for (iCycle = 0; iCycle < nCycles; iCycle++) {
      // Turn output LED ON
      digitalWrite(ledControlPin, HIGH);
      delay(ledOnTime);
      digitalWrite(ledControlPin, LOW);
      delay(ledOffTime);
    }

    // Turn them both off again
    digitalWrite(LED_BUILTIN, LOW);
    digitalWrite(ledControlPin, LOW);

    // reset gTrigger
    gTrigger = false;
  }
}

void triggerInterrupt() {
  if (~gTrigger) {
    // set the time that we received the trigger
    gTriggerMillis = millis();
    // set that we have turned it on
    gTrigger = true;
  }
}
