# $Id: Makefile,v 1.9 2008/11/20 20:01:07 jmuelmen Exp $

.PHONY: all
all:	help

.PHONY: help
help:	
	@echo "usage: make analysis NAME=<name of your new analysis>"

analysis:
	@if [ -z $(NAME) ] ; \
	then echo "usage: make analysis NAME=<name of your new analysis>" ;\
	else if [ -e $(NAME) ] ;\
		then echo "$(NAME) already exists" ;\
	else \
		mkdir $(NAME); \
		cp Templates/* $(NAME) ;\
		echo "Made a skeleton dilepton analysis in $(NAME)";\
		echo "";\
		echo "For a quick start, do the following:";\
		echo "   cd $(NAME)";\
		echo "   make -j 8 Results.tbl";\
		echo "";\
		echo "Cut names are defined in Looper.h.";\
		echo "Cuts are defined in Looper.cc.";\
		echo "Histograms are booked and filled in Looper.cc.";\
		echo "ConfigAndRun.cc defines samples and cuts to use for tables.";\
		echo "";\
		echo "To make a table, add a rule in ConfigAndRun.cc and run:";\
		echo "   make -j 8 <your table name>.tbl";\
	fi; fi