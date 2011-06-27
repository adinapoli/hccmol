#!/usr/bin/env python
# encoding: utf-8
"""
bioutils.py
Python bioinformatic scripting and everyday usage utilities.
Created by Alfredo Di Napoli on 2011-03-20.
"""

from Bio.PDB import *
from ftplib import FTP
from urllib2 import urlopen, URLError
from config import COMPOUNDS_DIR


ammino_tri2single = {
    'Ala' : 'A', 'Arg' : 'R', 'Asn' : 'N', 'Asp' : 'D', 'Cys' : 'C',
    'Gln' : 'Q', 'Glu' : 'E', 'Gly' : 'G', 'His' : 'H', 'Ile' : 'I',
    'Leu' : 'L', 'Lys' : 'K', 'Met' : 'M', 'Phe' : 'F', 'Pro' : 'P',
    'Ser' : 'S', 'Thr' : 'T', 'Trp' : 'W', 'Tyr' : 'Y', 'Val' : 'V'
}

ammino_single2tri = dict([[v,k] for k,v in ammino_tri2single.items()])


class SafePDBParser(PDBParser):
    """If the given file does not exists, it fetches it
    from the network."""
    
    def get_structure(self, name, filename):
        
        try:
            f = open("".join([COMPOUNDS_DIR,filename]), "r")
            
        except IOError:
            f = fetch_and_open(filename, "r")
            f.close()
            
        return PDBParser.get_structure(self, name, COMPOUNDS_DIR + filename)


def fetch_and_open(filename, mode):
    """Tries to open the file named filename. If it's not present into
    the local directory, attemps to fetch it from the network.
    Returns the opened stream"""
    
    f = None
    
    #The url where to get the amminoacid pdb file
    ammino_url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/pdb/"
    pdb_url = "http://www.rcsb.org/pdb/files/"
    
    #The compound name without extension
    compound_name = filename[:-4]
    
    
    try:
        f = open("".join([COMPOUNDS_DIR,filename]), mode)
    
    except IOError:
        
        query = compound_name.upper() + ".pdb"
        
        try:
            
            print "Downloading %s from PDBeChem repository..." % (filename)
            
            ftp = FTP('ftp.ebi.ac.uk')
            ftp.login()
            ftp.cwd("pub/databases/msd/pdbechem/files/pdb/")
            ftp.retrbinary('RETR ' + query, open(COMPOUNDS_DIR + query, 'w').write)
            ftp.close()
            
            f = open(COMPOUNDS_DIR + filename, mode)
            
            
        except Exception:
            
            try:
                
                print "Trying %s from PDB repository... " % (filename)
                
                web_url = urlopen(pdb_url + query)
                web_file = open(COMPOUNDS_DIR + filename, "w")
                web_file.write(web_url.read())
                web_file.close()
                web_url.close()
            
                f = open(COMPOUNDS_DIR + filename, mode)
                print ">>OK."
                
            except URLError:
                print query + " not found."
    
    return f
    
if __name__ == '__main__':
    f = fetch_and_open("ALA.pdb", "r")

