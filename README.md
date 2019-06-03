# EBI Test Data Generator
Computational approach to solving EBI, an algorithm used in computerized adaptive testing

## Use
* Allows user to generate test data with a single parameter varying
1) Follow menu prompt and select the appropriate option.
2) Enter the name of the output file. For example `theta_varies.csv`
3) Follow the remaining prompts for the selection

# Install needed libraries
* `python3 -m pip install --user numpy scipy matplotlib sympy`

# Caveat
In function `ebi()`, it is possible for `info_ts` to equal 0.
If `ebi() == 0` then a divide by 0 exception is raised.

A quick fix was to check if `info_ts` equals 0, and if it does then set it to a
valid number that is very very very close to 0, but allows the computation to succeed.
This was just a quick fix and a permanent fix will looked at in the future.
