#include <Keyboard.h>
#include <KeyboardLayout.h>

const int buttonPin = 4;         // input pin for pushbutton
int previousButtonState = HIGH;  // for checking the state of a pushButton
int counter = 0;                 // button push counter

void setup() {
  // make the pushButton pin an input:
  pinMode(buttonPin, INPUT);
  // initialize control over the keyboard:
  Keyboard.begin();
  // set LED to display
  pinMode(LED_BUILTIN, OUTPUT);
  // set previous button state
  previousButtonState = digitalRead(buttonPin);
  delay(100);
}

void loop() {
  // read the pushbutton:
  int buttonState = digitalRead(buttonPin);
  // if the button state has changed,
  if (buttonState != previousButtonState) {
      // and it's currently pressed:
    if (buttonState == HIGH) {
      // type out a message
      Keyboard.print("'");
      // turn LED on to show what is happening
      digitalWrite(LED_BUILTIN, HIGH);
      // delay, so as not to keep printing key
      delay(100);
    }
    else if (buttonState == LOW) {
      // turn LED off
      digitalWrite(LED_BUILTIN, LOW);
      delay(100);
    }
    // save the current button state for comparison next time:
    previousButtonState = buttonState;
  }
}
