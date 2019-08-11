#include <iostream>
#include <fstream>

#include "gcta.h"


ExprParserProb::ExprParserProb()
{
    symbol_table.add_variable("p",p);
    symbol_table.add_variable("q",q);
    
    //symbol_table.add_variable("y",y);
    //symbol_table.add_variable("z",z);
}


void ExprParserProb::parse_file(const string &filename)
{
    ifstream infile(filename);
    string line;
    unsigned count = 0;

    if (!infile.is_open())
        throw ("Error opening file " + filename + " for reading");

    cout << "##################" << endl;
    while ( getline (infile,line) )
    {
        if (line.size() == 0 || line[0] == '#') {
            continue;
        }

        cout << line << '\n';
        
        count++;
    }
    infile.close();
}