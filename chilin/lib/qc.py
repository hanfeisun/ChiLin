#!/usr/bin/env python
#Author     : meisl
#Name       :
#Description:
#Usage      :


class QC_Controller():
	"""
	All the class in the module derives from this class
	"""
	def __init__(self, template=""):
		self.template = template
		self.has_run = False
	
	def run(self):
		""" Run some QC tools or do some time-costing statistics """
		self.has_run = True
		return True

	def check(self):
		""" Check whether the quality of the dataset is ok. """
		if not self.has_run:
			self.run()
		return True

	def render(self, template = None):
		""" Generate the latex code for current section. """
		if not self.has_run:
			self.run()
		if template is None:
			template = self.template
		return ""
		
		


class RawQC(QC_Controller):
	def __init__(self):
		super(RawQC, self).__init__()
		print 'init basic_qc'
	def _basic_info(self):
		print 'basic QC information'
	def _fastqc_info(self):
		qc_main = ''
		files = []
		print 'fastqc\n'
	def run(self):
		self._basic_info()
		self._fastqc_info()
	def check(self):
		"""
		Check whether the FastQC's result is ok
		"""
		print 'if fastqc pass or not'
		
		
class MappingQC(QC_Controller):
	def __init__(self):
		super(MappingQC, self).__init__()		
		print 'init mapping qc'
	def _basic_mapping_statistics_info(self):
		print 'basic_mapping_statistics'
	def _mappable_ratio_info(self):
		print 'mappable_ratio'
	def _redundant_ratio_info(self):
		print 'redundant_ratio\n'
	def run(self):
		self._basic_mapping_statistics_info()
		self._mappable_ratio_info()
		self._redundant_ratio_info()
	def check():
		print 'mapping qc pass or not'


class PeakcallingQC(QC_Controller):
	def __init__(self):
		super(PeakcallingQC, self).__init__()
		print 'init peak calling  qc'
	def _peak_summary_info(self):
		print ' summary_info'

	def _velcro_ratio_info(self):
		print 'velcro_ratio_info '

	def _DHS_ratio_info(self):
		print 'DHS_ratio_info'

	def _replicate_info(self):
		print 'replicate_info\n'

	def run(self):
		self._peak_summary_info()
		self._velcro_ratio_info()
		self._DHS_ratio_info()
		self._replicate_info()
	def check():
		print 'pass or not'

		
class AnnotationQC(QC_Controller):
	def __init__(self):
		super(AnnotationQC, self).__init__()
		print 'intialization of function qc'
	def _ceas_info(self):
		print 'ceas qc'
	def _conservation_info(self):
		print 'conservation qc'
	def _motif_info(self):
		print 'motif info\n'
	def run(self):
		self._ceas_info()
		self._conservation_info()
		self._motif_info()
	def check(self):
		print 'pass or not'

