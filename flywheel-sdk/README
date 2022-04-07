Matlab SDK Installation Instructions

1. Unzip the flywheel-matlab-sdk zip file.
2. Add flywheel-sdk/api/rest-client.jar to the java classpath
    e.g. javaaddpath('flywheel-sdk/api/rest-client.jar')
3. Add flywheel-sdk to the matlab path
    e.g. addpath('flywheel-sdk')

Upgrade Note

When replacing the unzipped SDK, it's important to check that there is only
the most recent version of the rest-client.jar file on your javaclass path.

You can inspect your class path by running the command
  javaclasspath

This will produce a long list of java classes that Matlab includes, and at the end of the list you
should see a listing of the flywheel-sdk/api/rest-client.jar.
If you have more than one entry for the rest-client.jar, you can remove it using jmremove path.

For example:
  javarmpath('/full/path/to/outdated/rest-client.jar')
