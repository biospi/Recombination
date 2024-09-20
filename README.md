# Recombination finder
TODO
#https://github.com/YaolingYang/SparsePainter
https://github.com/SionBayliss/JGI_Seedcorn
https://github.com/danjlawson/finestructure4



sudo apt-get install cmake
sudo apt-get install libblas-dev liblapack-dev

/home/axel/tools/SparsePainter/SparsePainter 
-reffile /home/axel/python-projects/Recombination/output/wgs-mapping.tar.gz.filtered.ref.vcf 
-targetfile /home/axel/python-projects/Recombination/output/targets/5_SRR8204807.vcf 
-mapfile /home/axel/python-projects/Recombination/output/map.txt 
-popfile /home/axel/python-projects/Recombination/output/pop.txt 
-namefile /home/axel/python-projects/Recombination/output/targets/SRR8204807.txt
 -prob -haploid -chunklength -probstore raw