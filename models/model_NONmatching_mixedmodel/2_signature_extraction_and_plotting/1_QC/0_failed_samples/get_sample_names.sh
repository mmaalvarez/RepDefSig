
for sample in `grep FAI ../../1_parser_and_regressions/trace/trace.txt | cut -f2`
do
	tail -n 2 ../../1_parser_and_regressions/work/"$sample"*/.command.sh | head -n 1 | cut -d" " -f7 | sed "s/\[//g ; s/\]//g"
done > failed


for sample in `grep COM ../../1_parser_and_regressions/trace/trace.txt | cut -f2`
do
	tail -n 2 ../../1_parser_and_regressions/work/"$sample"*/.command.sh | head -n 1 | cut -d" " -f7 | sed "s/\[//g ; s/\]//g"
done > successful
