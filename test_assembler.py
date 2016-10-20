from unittest import TestCase
import Assemble

class TestAssembler(TestCase):

    def test_read_fasta(self):
        pass #self.fail()

    def test_example_works(self):
        """ Tests that small example works"""

        # Small given test case
        dbg_assembler = Assemble.Assembler(fragments_fasta=file("Data/coding_challenge_data_set_example.txt"))
        assembled_sequence = dbg_assembler.assemble()
        self.assertEqual(assembled_sequence, "ATTAGACCTGCCGGAATAC")


    def test_fail_on_euler_impossible(self):
        """ If it's impossible to make the graph Eulerian, output an understandable error"""

        non_eulerian_assembler = Assemble.Assembler(fragments_fasta=
                                                    file("Data/coding_challenge_data_set_example_made_non_eulerian.txt"))

        # Gives a reasonable error
        with self.assertRaises(AssertionError):
            non_eulerian_assembler.assemble()

