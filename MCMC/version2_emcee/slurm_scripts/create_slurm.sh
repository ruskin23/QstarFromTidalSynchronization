DIR=$1

for S in 10031409   10215422    10264202    10330495    10385682    10454401    10935310    10965961    1147276 11228612    11232745    11233911    11252617    11403216    11616200    11704044    12470530    2305543 3241344 3348093 3973504 4352168 4678171 4815612 4839180 5091614 5730394 7257373 7732791 7838639 7846730 8229048 8356054 8543278 8580438 8618226 8957954 10091257    10257903    10711551    10711913    1087921 10936427    10960995    10992733    11200773    11234677    11391181    12004679    12109845    12351927   2447893 2449090 2852560 2860788 3098194 3248019 3323289 3344427 3834364 3838496 4252226 4285087 4346875 4380283  4579321  4757331 4773155 4902030 4947726 5022440 5181455 5288543 5300878 5359678 5393558 5652260 5802470 5871918 6029130    6128027  6283224 6312521 6359798 6449552 6464285 6525196 6697716 6778050 6863840 6927629 6949550 6962018 7025540 7336754 7362852  7374746; 
do
./generate_slurm.sh $S $DIR
done
