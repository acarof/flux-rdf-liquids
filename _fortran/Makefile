all:compile

compile:
	gfortran -ffree-form -o gdr.exe gdr.f

compile_debug:
	gfortran -ffree-form -o gdr.exe gdr.f

run_debug:compile_debug
	./gdr.exe

regtest: clean run_debug
	diff -qr rdf_new.dat regtest/GDR_CM

clean:
	rm -rf *.exe rdf_new.dat
