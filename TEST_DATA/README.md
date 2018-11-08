Apply the vector fitting process to test data with

../bin/Vfit

Then follow the prompts for the data filename, model order, max iterations.

Alternatively use the script run_Vfit with command line options as follows:

run_Vfit filename model_order max_iternations

Examples:

1. Data from the impedance of a 'realistic' capacitor model consisting of a 
R L C series combination, R=85 mohms, L=43nH, C=1uF

run_Vfit CAPACITOR_CJS_Z.CSV 1 2

2. Measured impedance data for port 1 of a band pass filter with short circuit port 2.

run_Vfit BPF_SHORT_Z.CSV 12 4 
