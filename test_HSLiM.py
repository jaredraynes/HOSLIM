import pytest


@pytest.mark.sunsetimage
def test_HSLiM_output1():
	HSLiM_output =  h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw)

	data = [('as2-casein', 'B6VPY2', 'Domestic water buffalo', 10.3, 4, 3, 'FFIF', 6),
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 14.5, 8, 8, 'CLLAVALA', 15),
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 5.2, 3, 41, 'MAI', 43),
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 4.8, 3, 113, 'YLY', 115),
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 8.2, 3, 119, 'IVL', 121)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

def test_HSLiM_output2():
	HSLiM_output =  h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw)

	data = [('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 29.0, 14, 3, 'VLILACLVAVAVAM', 16),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 5.1, 3, 44, 'YYL', 46),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 8.1, 3, 69, 'LLL', 71),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 7.0, 3, 138, 'YFI', 140),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 5.4, 4, 143, 'AAVY', 146),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 6.6, 3, 151, 'LVY', 153),
			('CSN1 ECO0000313EMBLACU257801', 'D0QJ96', 'Duckbill platypus', 4.5, 3, 169, 'YAF', 171)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

	def test_HSLiM_output3():
	HSLiM_output =  h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw)

	data = [('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 28.5, 14, 3, 'VFIFACLVAVAMAV', 16),
			('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 6.8, 3, 33, 'LVM', 35),
			('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 6.1, 3, 43, 'ALV', 45),
			('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 5.4, 3, 62, 'MVY', 64),
			('D0QJA4_9MAMM', 'D0QJA4', 'Australian echidna', 9.5, 4, 74, 'YIFF',77)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

		def test_HSLiM_output4():
	HSLiM_output =  h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw)

	data = [('as2-casein', 'O97944', 'Dromedary', 9.9, 4, 3, 'FFIF', 6),
			('as2-casein', 'O97944', 'Dromedary', 15.7, 8, 8, 'CLLAVVLA', 15),
			('as2-casein', 'O97944', 'Dromedary', 6.0, 3, 42, 'VAI', 44),
			('as2-casein', 'O97944', 'Dromedary', 6.7, 3, 98, 'IVM', 100),
			('as2-casein', 'O97944', 'Dromedary', 6.5, 3, 167, 'FLW', 169)
			]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected

	#Also not working for the first HSLIM, the duplicator selector is not working properly, for some reason dropping if contains start = 3
	# e.g. for P02666 cow beta: gives 25.8	12	4	LILACLVALALA	15, but should be 28.6	13	3	VLILACLVALALA	15

#ADD TEST FOR P06796

#THIS WORKS?? WTF beta-casein	Q9GKK3	Horse	29.3	13	3	ILILACLVALALA	15
