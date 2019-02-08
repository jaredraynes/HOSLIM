import pytest


@pytest.mark.sunsetimage
def test_HSLiM_output1():
	HSLiM_output =  h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw)


	data = [('as2-casein', 'B6VPY2', 'Domestic water buffalo', 10.3, 4, 3, 'FFIF', 6),
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 14.5, 8, 8, 'CLLAVALA', 15)
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 5.2, 3, 41, 'MAI', 43)
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 4.8, 3, 113, 'YLY', 115)
			('as2-casein', 'B6VPY2', 'Domestic water buffalo', 8.2, 3, 119, 'IVL', 121)]
	df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
	expected = df.to_string(index=False)

	assert HSLiM_output == expected