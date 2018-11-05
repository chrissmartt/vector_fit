OBJS=read_data_from_file.o set_test_data.o update_poles.o calculate_residues.o calculate_fit_function.o calculate_error.o \
     write_data_to_plot.o write_coefficients.o eispack.o dmatrix.o

Vfit: Vfit.F90 $(OBJS)
	f95 -o Vfit Vfit.F90 $(OBJS)
        
read_data_from_file.o: read_data_from_file.F90
	f95 -c read_data_from_file.F90
        
set_test_data.o: set_test_data.F90
	f95 -c set_test_data.F90
        
update_poles.o: update_poles.F90
	f95 -c update_poles.F90
        
calculate_residues.o: calculate_residues.F90
	f95 -c calculate_residues.F90
                
calculate_fit_function.o: calculate_fit_function.F90
	f95 -c calculate_fit_function.F90
	
calculate_error.o: calculate_error.F90
	f95 -c calculate_error.F90
        
write_data_to_plot.o: write_data_to_plot.F90
	f95 -c write_data_to_plot.F90
        
write_coefficients.o: write_coefficients.F90
	f95 -c write_coefficients.F90
        
eispack.o: eispack.F90
	f95 -c eispack.F90
        
dmatrix.o: dmatrix.F90
	f95 -c dmatrix.F90

clean:
	rm Vfit $(OBJS)
