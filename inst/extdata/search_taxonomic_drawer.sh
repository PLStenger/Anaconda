#!/usr/bin/env line

# https://stackoverflow.com/questions/8101701/grep-characters-before-and-after-match?rq=1
# https://unix.stackexchange.com/questions/498027/grep-everything-up-until-and-including-a-pattern
# grep -oE '\S+Saksenaea' taxonomy.tsv
# https://unix.stackexchange.com/questions/24509/how-to-print-the-longest-line-in-a-file
# awk ' { if ( length > x ) { x = length; y = $0 } }END{ print y }'
# https://stackoverflow.com/questions/1521462/looping-through-the-content-of-a-file-in-bash


cat taxon_list.txt | while read line 
do
  #echo $line >> res.txt ;
   grep -oE '\S+'$line taxonomy.tsv | awk ' { if ( length > x ) { x = length; y = $0 } }END{ print y }' >> res.txt
done

# Delete duplicates lines
awk '!seen[$0]++' res.txt > res_02.txt
sed -E '/^$/d' res_02.txt > res_03.txt
sed 's/taxonomy.tsv://g' res_03.txt > taxon_list_drawer.txt


rm res.txt
rm res_02.txt
rm res_03.txt