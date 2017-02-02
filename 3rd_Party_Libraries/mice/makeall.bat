rem This script builds the SPICE delivery
rem for the mice package of the toolkit.
rem
rem The script must be executed from the
rem mice directory.
rem
cd src
rem
rem Creating cspice
rem
cd cspice
call mkprodct.bat
cd ..
rem
rem Creating csupport
rem
cd csupport
call mkprodct.bat
cd ..
rem
rem Creating mice
rem
cd mice
call mkprodct.bat
cd ..
rem
rem Creating brief_c
rem
cd brief_c
call mkprodct.bat
cd ..
rem
rem Creating chrnos_c
rem
cd chrnos_c
call mkprodct.bat
cd ..
rem
rem Creating ckbref_c
rem
cd ckbref_c
call mkprodct.bat
cd ..
rem
rem Creating commnt_c
rem
cd commnt_c
call mkprodct.bat
cd ..
rem
rem Creating cook_c
rem
cd cook_c
call mkprodct.bat
cd ..
rem
rem Creating frmdif_c
rem
cd frmdif_c
call mkprodct.bat
cd ..
rem
rem Creating inspkt_c
rem
cd inspkt_c
call mkprodct.bat
cd ..
rem
rem Creating mkspk_c
rem
cd mkspk_c
call mkprodct.bat
cd ..
rem
rem Creating msopck_c
rem
cd msopck_c
call mkprodct.bat
cd ..
rem
rem Creating spacit_c
rem
cd spacit_c
call mkprodct.bat
cd ..
rem
rem Creating spkdif_c
rem
cd spkdif_c
call mkprodct.bat
cd ..
rem
rem Creating spkmrg_c
rem
cd spkmrg_c
call mkprodct.bat
cd ..
rem
rem Creating tobin_c
rem
cd tobin_c
call mkprodct.bat
cd ..
rem
rem Creating toxfr_c
rem
cd toxfr_c
call mkprodct.bat
cd ..
rem
rem Creating versn_c
rem
cd versn_c
call mkprodct.bat
cd ..
rem
rem Creating micecook
rem
cd micecook
cd ..
cd ..
rem Toolkit Build Complete
