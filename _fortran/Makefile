WFLAGS   = -Waliasing -Wampersand -Wc-binding-type -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-tabs -Wtarget-lifetime -Wunderflow -Wunused-but-set-variable -Wunused-variable -Wconversion -Werror


all:compile

compile:
	gfortran -ffree-form -o gdr.exe gdr.f

compile_debug:
	gfortran -ffree-form  $(WFLAGS) -o gdr.exe gdr.f

run_debug:compile_debug
	./gdr.exe

regtest: clean run_debug
	diff -qr rdf_new.dat regtest/rdf_new_d7a67cc.dat

clean:
	rm -rf *.exe rdf_new.dat
