perl _format_blast_GO_db.pl > _grassesAthaGOdb.faa

makeblastdb -in _grassesAthaGOdb.faa -dbtype prot 

rm _grassesAthaGOdb.faa
mv _grassesAthaGOdb.faa.p* ..
