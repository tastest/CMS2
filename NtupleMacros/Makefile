# $Id: Makefile,v 1.12 2008/12/05 16:28:49 jmuelmen Exp $

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
		mkdir -p $(NAME); \
		find Templates/ -maxdepth 1 -type f -exec sh -c "cp {} $(NAME)" \; ;\
		echo "Made a skeleton dilepton analysis in $(NAME)" | tee $(NAME)/AAAREADME.txt ;\
		echo "" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "For a quick start, do the following:" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "   cd $(NAME)" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "   make -j 8 Results.tbl" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "Cut names are defined in Looper.h." | tee -a $(NAME)/AAAREADME.txt ;\
		echo "Cuts are defined in Looper.cc." | tee -a $(NAME)/AAAREADME.txt ;\
		echo "Histograms are booked and filled in Looper.cc." | tee -a $(NAME)/AAAREADME.txt ;\
		echo "ConfigAndRun.h defines samples and cuts to use for tables." | tee -a $(NAME)/AAAREADME.txt ;\
		echo "" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "To make a table, add a rule in ConfigAndRun.h and run:" | tee -a $(NAME)/AAAREADME.txt ;\
		echo "   make -j 8 <your table name>.tbl" | tee -a $(NAME)/AAAREADME.txt ;\
	fi; fi
