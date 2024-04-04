#!/bin/sh

#outfile=output.txt
#outfile=output_1000bins.txt
#outfile=output2.txt
#rm -f $outfile
#make
#./cadmix_mt --mutations 500 -b 50 -t 0 -o
#./cadmix_mt --mutations 51
#./cadmix_mt --mutations 500 -b 100 -t 0 -o op.s500.b100.txt
#./cadmix_mt -a 0.02 -m 3500 -b 2

#./cadmix_mt -a 0.1 -m 500 -b 40 -t 8 -o op.s500.b40.a01.test2.txt

#./cadmix_mt -a 0.1 -m 500 -b 18 -t 4 -o op.s500.b18.a01.test2.txt
#./cadmix_mt -a 0.02 --fast -m 500 -b 20 -t 4
#./cadmix_mt -a 0.1 --fast  -p 2000 -m 5000 -b 20 -t 4
#./cadmix_mt -a 0.02 -p 2000 -m 100000 -b 10 -A 4 -B 4 --fast 

#./cadmix_mt -a 0.02 -p 2000 -m 50 -w 100000 -A 2 --fast
#./cadmix_mt -a 0.02 -p 2000 -m 100000 -A 2 -b 20 -t 4 --fast
#./cadmix_mt -a 0.02 -p 2000 -m 100000 -A 2 -b 20 -t 4 --fast -B 3


#./cadmix_mt -a 0.02 -p 2000  -m 10000 -b 200 -t 8 -A 1000 -B 12000 --cond --non-ind "neanderthal_nonind_p2k_corret?.txt"
#-o "LengthCond/neanderthal_nonind_p10k.txt"
#./cadmix_mt -a 0.02 -p 1000  -m 10000 -b 200 -t 8 -A 1000 -B 12000 --cond --non-ind  --verbose
#-o "LengthCond/neanderthal_ind_p10k.mult.txt" --mult-ref


#./cadmix_mt -a 0.8 -p 12000 -P 8000 -m 10000 -w 20000 -A 10 -B 2000 --mult-ref --cond -u $u -n $n | awk 'NR==8 || NR==20'
#outf=newneanderthalrun_first.txt
#rm -f ${outf}
#./cadmix_mt -t 4 -b 20 -a 0.01 -p 2000 -m 100000 -A 500 -B 2000 -c 10 --cond --out ${outf}

#./cadmix_mt -a 0.8 -p 12000 -P 8000 -m 10000 -b 200 -t 8 -A 10 -B 2000 --mult-ref --cond 
#for w in 2000, 20000, 100000, 1000000, 10000000
#do
#    for n in 1, 10, 100
#    do    
#        echo $w $n
#        ./cadmix_mt -a 0.8 -p 12000 -P 8000 -m 10000 -w $w -A 10 -B 2000 --mult-ref -n $n | tail -n 5 | head -n 2
#    done
#done
./cadmix_mt -a 0.1 -p 8000 -P 12000 -m 100 -A 10 -B 2000 -e 10000 --cond --w 20000 -u 1.5E-9
#./cadmix_mt -a 0.8 -p 12000 -P 8000 -m 100000 -w 10000 -A 10 -B 2000 --mult-ref --cond -n 100

#./cadmix_mt -a 0.02 -p 2000 -w 100000 -A 2 --fast

#./cadmix_mt -a 0.02 -p 2000 -m 100000 -b 40 -t 4 -B 4 -A 4 --fast -o op.s100k.b80.a002.p2k.fast.txt

#./cadmix_mt -a 0.5 -m 100000 -B 3 -A 4 -b 100 -t 4 --fast 
#-o 
#op.s100k.a05.b100.tb2000.ta10.txt

#-o op.s1000.b20.a002.txt
#./cadmix_mt -b 20 -t -o op.s500.b20.exact.txt
#./cadmix_mt
#./cadmix_mt --mutations 500 -b 20 -t -o op.s500.b20.exact.txt
# -t -p op.s500.b20.txt
#-b 2 -o op.extra.1.txt

#-o output_exact.200.txt

#outfile=op.s100k.b100.a002.fast.txt
#outfile=./outputs/op.s100k.b100.a002.2.txt
#outfile2=./outputs/op.s100k.b100.a002.fast.txt
#outfile=op.s100k.a05.b100.tb2000.ta10.txt

outfile=op.s100k.b100.a002.p2k.txt
outfile2=op.s100k.b100.a002.p2k.fast.txt

#./makeplot.py -i $outfile -i2 $outfile2 -nt1 4 -nr 10
