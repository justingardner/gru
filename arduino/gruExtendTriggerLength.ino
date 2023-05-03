/*
  Function to extend the trigger length of digital pulse from magnet so that
  we can read it from our NI board
 */

volatile int state = LOW;
int outputPin = 13;
// length of output trigger in milliseconds
int gOutTriggerLen = 45;
// time of last trigger
unsigned long gTriggerMillis = 0;
// triggered
int gTrigger = false;

void setup() {
  // initialize digital pin 13 as an output.
  pinMode(outputPin, OUTPUT);
  // attach the interrupt to detect digital pulses
  // to interrupt 0 which is line 2 on the Uno
  attachInterrupt(0, triggerIn, RISING);
  // set the initial output to low
  digitalWrite(outputPin,LOW);
}

// the loop function runs over and over again forever
void loop() {
  // check if we are past the time needed
  // to trun the output trigger off. Note that
  // the cast to unsigned long should make this handle
  // rollover of the millis function correctly
  if ((gTrigger) && ((unsigned long)(millis()-gTriggerMillis)>=gOutTriggerLen)) {
    // turn off output pin
    digitalWrite(outputPin, LOW); 
    // and set that we are no longer triggered
    gTrigger = false;
  }
}

void triggerIn()
{
  if (~gTrigger) {
    // set the time that we received the trigger
    gTriggerMillis = millis();
    // set output pin high
    digitalWrite(outputPin, HIGH);
    // set that we have turned it on
    gTrigger = true;
  }
}
