/***************** start of the getcgivars() module **********************/

/*************************************************************************/
/**                                                                     **/
/**     getcgivars.c-- routine to read CGI input variables into an      **/
/**         array of strings.                                           **/
/**                                                                     **/
/**     Written in 1996 by James Marshall, james@jmarshall.com, except  **/
/**     that the x2c() and unescape_url() routines were lifted directly **/
/**     from NCSA's sample program util.c, packaged with their HTTPD.   **/
/**                                                                     **/
/**     For the latest, see http://www.jmarshall.com/easy/cgi/ .        **/
/**                                                                     **/
/*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/** Convert a two-char hex string into the char it represents. **/
char x2c(char *what) {
  register char digit;

  digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
  return(digit);
}

/** Reduce any %xx escape sequences to the characters they represent. **/
void unescape_url(char *url) {
  register int i,j;

  for(i=0,j=0; url[j]; ++i,++j) {
    if((url[i] = url[j]) == '%') {
      url[i] = x2c(&url[j+1]) ;
      j+= 2 ;
    }
  }
  url[i] = '\0' ;
}


/** Read the CGI input and place all name/val pairs into list.        **/
/** Returns list containing name1, value1, name2, value2, ... , NULL  **/
char **getcgivars() {
  register int i ;
  char *request_method ;
  int content_length;
  char *cgiinput ;
  char **cgivars ;
  char **pairlist ;
  int paircount ;
  char *nvpair ;
  char *eqpos ;

  /** Depending on the request method, read all CGI input into cgiinput. **/
  request_method= getenv("REQUEST_METHOD") ;

  if (!strcmp(request_method, "GET") || !strcmp(request_method, "HEAD") ) {
    /* Some servers apparently don't provide QUERY_STRING if it's empty, */
    /*   so avoid strdup()'ing a NULL pointer here.                      */
    char *qs ;
    qs= getenv("QUERY_STRING") ;
    cgiinput= strdup(qs  ? qs  : "") ;
  }
  else if (!strcmp(request_method, "POST")) {
    /* strcasecmp() is not supported in Windows-- use strcmpi() instead */
    if ( strcasecmp(getenv("CONTENT_TYPE"), "application/x-www-form-urlencoded")) {
      printf("Content-Type: text/plain\n\n") ;
      printf("getcgivars(): Unsupported Content-Type.\n") ;
      exit(1) ;
    }
    if ( !(content_length = atoi(getenv("CONTENT_LENGTH"))) ) {
      printf("Content-Type: text/plain\n\n") ;
      printf("getcgivars(): No Content-Length was sent with the POST request.\n") ;
      exit(1) ;
    }
    if ( !(cgiinput= (char *) malloc(content_length+1)) ) {
      printf("Content-Type: text/plain\n\n") ;
      printf("getcgivars(): Couldn't malloc for cgiinput.\n") ;
      exit(1) ;
    }
    if (!fread(cgiinput, content_length, 1, stdin)) {
      printf("Content-Type: text/plain\n\n") ;
      printf("getcgivars(): Couldn't read CGI input from STDIN.\n") ;
      exit(1) ;
    }
    cgiinput[content_length]='\0' ;
  }
  else {
    printf("Content-Type: text/plain\n\n") ;
    printf("getcgivars(): Unsupported REQUEST_METHOD.\n") ;
    exit(1) ;
  }

  /** Change all plusses back to spaces. **/
  for (i=0; cgiinput[i]; i++) if (cgiinput[i] == '+') cgiinput[i] = ' ' ;

  /** First, split on "&" and ";" to extract the name-value pairs into **/
  /**   pairlist.                                                      **/
  pairlist= (char **) malloc(256*sizeof(char **)) ;
  paircount= 0 ;
  nvpair= strtok(cgiinput, "&;") ;
  while (nvpair) {
    pairlist[paircount++]= strdup(nvpair) ;
    if (!(paircount%256))
      pairlist= (char **) realloc(pairlist,(paircount+256)*sizeof(char **)) ;
    nvpair= strtok(NULL, "&;") ;
  }
  pairlist[paircount]= 0 ;    /* terminate the list with NULL */

  /** Then, from the list of pairs, extract the names and values. **/
  cgivars= (char **) malloc((paircount*2+1)*sizeof(char **)) ;
  for (i= 0; i<paircount; i++) {
    if (eqpos=strchr(pairlist[i], '=')) {
      *eqpos= '\0' ;
      unescape_url(cgivars[i*2+1]= strdup(eqpos+1)) ;
    } else {
      unescape_url(cgivars[i*2+1]= strdup("")) ;
    }
    unescape_url(cgivars[i*2]= strdup(pairlist[i])) ;
  }
  cgivars[paircount*2]= 0 ;   /* terminate the list with NULL */

  /** Free anything that needs to be freed. **/
  free(cgiinput) ;
  for (i=0; pairlist[i]; i++) free(pairlist[i]) ;
  free(pairlist) ;

  /** Return the list of name-value strings. **/
  return cgivars ;

}

/***************** end of the getcgivars() module ********************/

#include "../xtalcomp.h"

#include <iostream>
#include <sstream>
#include <string>

#define PRINT_DIV \
  printf("|-------------------------------------|-------------------------------------|\n")

bool parsePOSCAR(char *, XcMatrix&, std::vector<XcVector>&, std::vector<unsigned int>&);

std::string debug;

// XtalComp CGI wrapper:
int main() {
  char ** cgivars = getcgivars();

  XcMatrix cell1 (4.5, 0.0, 0.0,
                  1.2, 2.4, 0.0,
                  2.5, 6.4, 1.1);

  XcMatrix cell2 (2.5, 0.0, 0.0,
                  2.3, 4.2, 0.0,
                  1.5, 2.1, 2.1);

  std::vector<XcVector> pos1 (5);
  std::vector<XcVector> pos2 (5);

  std::vector<unsigned int> types1 (5);
  std::vector<unsigned int> types2 (5);

  float cartTol, angleTol;

  float transform[16];

  bool validInput = true;

  for (int i=0; cgivars[i]; i+= 2) {
    char *key = cgivars[i];
    char *val = cgivars[i+1];
    if (strcmp(key, "pos1") == 0) {
      if (!parsePOSCAR(val, cell1, pos1, types1)) validInput = false;
    }
    else if (strcmp(key, "pos2") == 0) {
      if (!parsePOSCAR(val, cell2, pos2, types2)) validInput = false;
    }
    else if (strcmp(key, "cartTol") == 0) {
      if (sscanf(val, "%f", &cartTol) != 1) validInput = false;
    }
    else if (strcmp(key, "angleTol") == 0) {
      if (sscanf(val, "%f", &angleTol) != 1) validInput = false;
    }
  }

  if (!validInput) {
    printf("Content-type: text/html\n\n");
    printf("<html>\n");
    printf("<head><title>XtalComp Results</title></head>\n");
    printf("<body>\n");
    printf("<h1>Invalid input</h1>\n");
    printf("Go back and check your inputs.\n");
    printf("<br>");
    //printf(debug.c_str());
    printf("</body>\n");
    printf("</html>\n");
    return 1;
  }

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 transform, cartTol, angleTol);

  printf("Content-type: text/html\n\n");

  printf("<html>\n");
  printf("<head><title>XtalComp Results</title></head>\n");
  printf("<body>\n");
  printf("<h1>Result:</h1>\n");
  printf("Using a cartesian tolerance of %f and an angular tolerance of %f...<br><br>\n",
         cartTol, angleTol);
  printf("The structures %s match!<br>\n", (match) ? "DO" : "do NOT");
  if (match) { // Print transform
    printf("<font face=\"Courier New, Courier, monospace\">\n");
    printf("<pre>\n");
    printf("<h2>Transformation matrix:</h2>\n");
    printf("|--%10s--%10s--%10s--%10s--|\n", "----------", "----------",
           "----------", "----------");
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[0*4+0], transform[0*4+1], transform[0*4+2], transform[0*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[1*4+0], transform[1*4+1], transform[1*4+2], transform[1*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[2*4+0], transform[2*4+1], transform[2*4+2], transform[2*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[3*4+0], transform[3*4+1], transform[3*4+2], transform[3*4+3]);
    printf("|--%10s--%10s--%10s--%10s--|\n", "----------", "----------",
           "----------", "----------");
    printf("</pre>\n");
    printf("</font>\n") ;
    }

  printf("<h1>Input structures:</h1>\n") ;
  printf("<font face=\"Courier New, Courier, monospace\">\n");
  printf("<pre>\n");
  PRINT_DIV;
  printf("| %-35s | %-35s |\n",
         "First cell matrix (row vectors)",
         "Second cell matrix (row vectors)");
  PRINT_DIV;
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[0][0], cell1[0][1], cell1[0][2], "",
         cell2[0][0], cell2[0][1], cell2[0][2], "");
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[1][0], cell1[1][1], cell1[1][2], "",
         cell2[1][0], cell2[1][1], cell2[1][2], "");
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[2][0], cell1[2][1], cell1[2][2], "",
         cell2[2][0], cell2[2][1], cell2[2][2], "");
  PRINT_DIV;
  printf("| %-35s | %-35s |\n",
         "type: fractional coordinate",
         "type: fractional coordinate");
  PRINT_DIV;
  for (int i = 0; i < types1.size(); ++i) {
    printf("| %3hu: %9.5f %9.5f %9.5f %0s | %3hu: %9.5f %9.5f %9.5f %0s |\n",
           types1[i], pos1[i][0], pos1[i][1], pos1[i][2], "",
           types2[i], pos2[i][0], pos2[i][1], pos2[i][2], "");
  }
  PRINT_DIV;

  //printf("debug: \n%s\n", debug.c_str());
  printf("</pre>\n");
  printf("</font>\n") ;
  //printf("<!--#include virtual=\"xtalcompcounter.cgi\" -->\n");
  printf("Comparisons performed since June 6th, 2011: "
         "<br><embed src=\"xtalcompcounter.cgi\" />\n");
  printf("</body>\n") ;
  printf("</html>\n") ;

  /** Free anything that needs to be freed **/
  for (int i=0; cgivars[i]; i++) free(cgivars[i]) ;
  free(cgivars) ;

  exit(0) ;
}

void Debug(const char *str, const double d)
{
  char buffer[128];
  snprintf(buffer, 32, "%s %f\n", str, d);
  debug += buffer;
}
void Debug(const std::string &str, const double d) {Debug(str.c_str(), d);}

bool parsePOSCAR(char *str, XcMatrix &cell,
                 std::vector<XcVector> &pos,
                 std::vector<unsigned int> &types)
{
  std::string stdstr (str);
  std::istringstream lines (stdstr);
  std::string line;
  bool cart = false;

  // First line is comment
  getline(lines, line);

  // Next line is scale factor
  getline(lines, line);
  float scale;
  if (sscanf(line.c_str(), "%f", &scale) != 1) return false;

  // Next comes the matrix
  float x,y,z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[0][0] = x;
  cell[0][1] = y;
  cell[0][2] = z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[1][0] = x;
  cell[1][1] = y;
  cell[1][2] = z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[2][0] = x;
  cell[2][1] = y;
  cell[2][2] = z;

  // Apply scale:
  cell *= scale;

  // Store frac->cart matrix
  XcMatrix toCart = cell.transpose().inverse();

  // List of atom types
  std::vector<int> counts (15); // Allow up to 15 atom types.
  getline(lines, line);
  int tmpint;
  int numTypes = sscanf(line.c_str(), "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                        &counts[0], &counts[1], &counts[2],
                        &counts[3], &counts[4], &counts[5],
                        &counts[6], &counts[7], &counts[8],
                        &counts[9], &counts[10], &counts[11],
                        &counts[12], &counts[13], &counts[14], &tmpint);

  if (numTypes > 15) return false;

  // Starts with either [Ss]elective dynamics, [KkCc]artesian, or
  // other for fractional coords.
  getline(lines, line);

  // If selective dynamics, get the next line
  if (line.at(0) == 'S' || line.at(0) == 's')
    getline(lines, line);

  // Check if we're using cartesian or fractional coordinates:
  if (line.at(0) == 'K' || line.at(0) == 'k' ||
      line.at(0) == 'C' || line.at(0) == 'c' )
    cart = true;
  else
    cart = false;


  // Coordinates
  // determine number of atoms:
  types.clear();
  int numAtoms = 0;
  for (int i = 0; i < numTypes; ++i) {
    numAtoms += counts[i];
    for (int j = 0; j < counts[i]; ++j) {
      types.push_back(i);
    }
  }

  types.resize(numAtoms);

  Debug("numAtoms:", numAtoms);

  // Grab vectors
  XcVector tmp;
  pos.clear();
  for (int atom_i = 0; atom_i < numAtoms; ++atom_i) {
    getline(lines, line);
    if (sscanf(line.c_str(), "%f %f %f",
               &x, &y, &z) != 3) return false;
    tmp = XcVector(x,y,z);
    debug += "pos line: " + line + "\n";
    Debug("x: ", tmp.x());
    Debug("y: ", tmp.y());
    Debug("z: ", tmp.z());
    if (cart) {
      tmp = toCart * tmp;
      debug += "Converted to cartesian:\n";
      Debug("x: ", tmp.x());
      Debug("y: ", tmp.y());
      Debug("z: ", tmp.z());
    }
    pos.push_back(tmp);
  }

  Debug("pos size: ", pos.size());
  Debug("types size: ", types.size());
  return true;
}
