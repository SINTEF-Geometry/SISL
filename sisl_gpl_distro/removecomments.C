//===========================================================================
//                                                                           
// File: removecomments.C                                                    
//                                                                           
// Created: Thu Feb 10 13:45:40 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: removecomments.C,v 1.1 2005-02-25 15:02:08 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <iostream>

using namespace std;

int main()
{
    char c;
    bool in_comment = false;
    bool info_on_line = false;
    bool comment_in_this_line = false;
    do {
	cin.read(&c,1);
	if (in_comment) {
	    // looking for end of comment
	    int asterix_count = 0;
	    while (c == '*' && !cin.eof()) {
		++asterix_count;
		cin.read(&c, 1);
	    }
	    if (asterix_count > 0 && c == '/') {
		// this is the end of the comment
		in_comment = false;
	    }
	} else {
	    int slash_count = 0;
	    while (c == '/' && !cin.eof()) {
		++slash_count;
		cin.read(&c, 1);
	    }
	    
	    if (slash_count > 0 && c == '*') {
		// new comment starts here
		in_comment = true;
		comment_in_this_line = true;
		--slash_count; // the last slash should not be written
	    }
	    for (int i = 0; i < slash_count; ++i) {
		cout << '/';
	    }
	    if (!in_comment) {
		if (c == '\n') {
		    if (info_on_line || !comment_in_this_line) {
			cout.write(&c, 1);
		    }
		    info_on_line = false;
		    comment_in_this_line = false;
		} else {
		    info_on_line = true;
		    cout.write(&c, 1);
		}
	    }
	}
    } while (!cin.eof());
    //while (c != EOF);

    return 0;
};
