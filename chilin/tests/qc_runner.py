from chilin.qc import (
	RawQC,
	MappingQC,
	PeakcallingQC,
	AnnotationQC)
import unittest




class RawQCTestCase(unittest.TestCase):
    def setUp(self):
        self.rawqc = RawQC()
    def testRun(self):
        self.rawqc.run()



if __name__ == '__main__':
	unittest.main()

		
	
	
