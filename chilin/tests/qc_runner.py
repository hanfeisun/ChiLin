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

class MappingQCTestCase(unittest.TestCase):
    def setUp(self):
        self.mappingqc = MappingQC()
    def testRun(self):
        self.mappingqc.run()
        
class PeakcallingQCTestCase(unittest.TestCase):
    def setUp(self):
        self.peakcallingqc = PeakcallingQC()
    def testRun(self):
        self.peakcallingqc.run()

class AnnotationQCTestCase(unittest.TestCase):
    def setUp(self):
        self.annotationqc = AnnotationQC()
    def testRun(self):
        self.annotationqc.run()


if __name__ == '__main__':
	unittest.main()

		
	
	
