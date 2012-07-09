#!/usr/bin/env python
#Author     : meisl
#Name       :
#Description:
#Usage      :
#------------------------------

#------------------------------
class QC():
	def __init__(self):
		print 'df'

class Basic_qc(QC):
	def __init__(self):
		print 'init basic_qc'
	def __basic_info(self):
		print 'basic QC information'
	def __fastqc_info(self):
		qc_main = ''
		files = []
		print 'fastqc\n'
	def summary(self):
		self.__basic_info()
		self.__fastqc_info()
	def check(self):
		print 'if fastqc pass or not'	
#------------------------------
class Mapping_qc(QC):
	def __init__(self):
		print 'init mapping qc'
	def __basic_mapping_statistics_info(self):
		print 'basic_mapping_statistics'
	def __mappable_ratio_info(self):
		print 'mappable_ratio'
	def __redundant_ratio_info(self):
		print 'redundant_ratio\n'
	def summary(self):
		self.__basic_mapping_statistics_info()
		self.__mappable_ratio_info()
		self.__redundant_ratio_info()
	def check():
		print 'mapping qc pass or not'

#-------------------------------
class Peak_calling_qc(QC):
	def __init__(self):
		print 'init peak calling  qc'
	def __peak_summary_info(self):
		print ' summary_info'

	def __velcro_ratio_info(self):
		print 'velcro_ratio_info '

	def __DHS_ratio_info(self):
		print 'DHS_ratio_info'

	def __replicate_info(self):
		print 'replicate_info\n'

	def summary(self):
		self.__peak_summary_info()
		self.__velcro_ratio_info()
		self.__DHS_ratio_info()
		self.__replicate_info()
	def check():
		print 'pass or not'
#------------------------------
class Function_qc(QC):
	def __init__(self):
		print 'nitialization of function qc'
	def __ceas_info(self):
		print 'ceas qc'
	def __conservation_info(self):
		print 'conservation qc'
	def __motif_info(self):
		print 'motif info\n'
	def summary(self):
		self.__ceas_info()
		self.__conservation_info()
		self.__motif_info()
	def check(self):
		print 'pass or not'

