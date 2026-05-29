cp ./7.rp.qmd ./rp_final
cd rp_final
pixi run quarto render ./7.rp.qmd 
rm ./7.rp.qmd 
tar -czf ../16S_suran.tar.gz ./*
zip -r ../16S_suran.zip ./* -x "*.DS_Store"