
function [s, fSerial] = triggerBoxTest(s, fSerial)

if nargin == 0
  % display all serial devices
  serialDev = dir('/dev/cu.*');

  nSerialDev = length(serialDev);
  for i = 1:nSerialDev
    disp(sprintf('%i: %s', i,serialDev(i).name));
  end
  serialDevNum = getnum('Choose a serial device (0 to quit)',0:length(serialDev));
  if serialDevNum == 0,return,end

  % try to open the device
  try
    s = serial(fullfile('/dev',serialDev(serialDevNum).name));
  
  baudRate = 57600;
  set(s,'BaudRate',baudRate);
  %set(s,'Parity',parity);
  %set(s,'StopBits',stopBits);
  %set(s,'DataBits',dataLen);
  %set(s,'Terminator','CR/LF');
  %set(s,'InputBufferSize',2048);
    fopen(s);
    fSerial = s;
  catch
    disp(sprintf('(moncalib:initSeriaPortUsingSerial) Failed to open serial port. Sometimes restarting matlab helps. Or plugging in and unplugging in serial device. Or making sure that you have the serial device driver installed.'));
  end
  % get help
  fwrite(s,'[?]');
  str = '';
  while (s.BytesAvailable ~= 0)
    str = strcat(str,fread(s,s.BytesAvailable)');
  end
  disp(sprintf('%s',str));
end
