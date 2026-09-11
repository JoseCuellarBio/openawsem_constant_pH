from openawsem.helperFunctions.generate_charge import generate_charge_file


def test_generate_charge_file_uses_constant_ph_defaults(tmp_path):
    fasta = tmp_path / "protein.fasta"
    fasta.write_text(">test\nDECKRHYAG\n")

    output = tmp_path / "charge_pH.txt"
    generate_charge_file(fasta, output)

    assert output.read_text().splitlines() == [
        "0 -1.0",
        "1 -1.0",
        "2 -1.0",
        "3 1.0",
        "4 1.0",
        "5 1.0",
        "6 -1.0",
        "7 0.0",
        "8 0.0",
    ]
