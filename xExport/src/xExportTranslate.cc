/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
*/

#include "xExportTranslate.h"
#define TAILLE_CHARAC_ILLEGAL 7

using namespace std;
namespace xexport
{
xExportTranslate::xExportTranslate(string b, string src, string nom)
    : xExportEnsightBinary(), project_name_(nom), cSrc_(src), isBinary_(isFileBinary(b))
{
   if (isBinary_)
      fpos_.open(src.c_str(), ifstream::in);
   else
      fpos_.open(src.c_str(), ifstream::in | ifstream::binary);
   if (!fpos_.is_open())
   {
      std::cerr << "Erreur : Ouverture impossible, le nom ou le chemin n'est pas correct" << std::endl;
      exit(1);
   }
   readFileGmsh();
}
bool xExportTranslate::isFileBinary(const string& test) { return (test[0] == 'b') || (test[0] == 'B'); }

void xExportTranslate::readFileGmsh()
{
   if (isBinary_)
   {
      readGmshBinary();
   }
   else
   {
      readGmshAscii();
   }

   writeFileEnsight();
}

void xExportTranslate::readGmshAscii()
{
   string ligne;
   // float x1,y,z;
   // CommenÃ§ons Ã  lire le fichier ligne par ligne et a rentrer les information.

   getline(fpos_, ligne);
   // obtention de la premiÃšre ligne

   // fig.addNomFigure(giveNomProject(ligne));
   fig.addNomFigure(project_name_);
   // rÃ©cupÃ©ration du type d'Ã©lÃ©ment

   getline(fpos_, ligne);

   // le j c'est la rÃ©fÃ©rence

   int j = ligne.find_first_of("(");

   string test = ligne.substr(0, j);
   fig.addTypeElem(giveTypeElementGmsh(test[1]));
   fig.addtypeVariable(giveTypeVariable(test[0]));
   // attention il y a des types complexe que je ne traite pas pour le moment.
   string x("");
   // sert a la conversion
   std::istringstream ins;
   do
   {
      int b = ligne.find_first_of(")");
      j = ligne.find_first_of("(") + 1;
      float hj[3];
      float bte;
      int cpt = -1;
      for (; j <= b; j++)
      {
         if (!isSeparateur(ligne[j]) && !isFinAc(ligne[j]))
         {
            x += ligne[j];
         }
         else
         {
            cpt++;

            ins.str(x);
            ins >> bte;
            hj[cpt] = bte;
            if (cpt == 2)
            {
               fig.addElem(hj[0], hj[1], hj[2]);

               cpt = -1;
            }

            ins.clear();
            x.erase();
         }
      }  // Debut travail sur Variables

      traiteVariable(ligne.substr(ligne.find_first_of('{')));

   } while (getline(fpos_, ligne) && ligne[0] != '}');

   fpos_.close();
}

short xExportTranslate::giveTypeElementGmsh(char a)
{
   switch (a)
   {
      case 'S':
         return 4;
         break;  // quadrangle4
      case 'T':
         return 3;
         break;  // triangle 3
      case 'P':
         return 1;
         break;  // point
      case 'L':
         return 2;
         break;  // line
      case 'Q':
         return 5;
         break;  // tetra
      case 'H':
         return 8;
         break;  // hexa
         // on continue comme pour faire la conversion
      default:
         std::cerr << "erreur de type" << std::endl;
         exit(1);
         break;
   }
   // si return 0, erreur;
}

short xExportTranslate::giveTypeVariable(const char nan)
{
   switch (nan)
   {
      case 'S':
         return 1;
         break;
      case 'V':
         return 2;
         break;
      case 'T':
         return 3;
         break;
      default:
         throw;
         return 0;
   }
}

void xExportTranslate::writeFileEnsight()
{
   openFile(fig.getProjectName());
   startView(fig.getProjectName());
   endView();
   closeFile();
}

/******************************************************************/

string xExportTranslate::giveNomProject(string li)
{
   int premiere;
   int derniere;
   size_t i, j;
   char illegal[TAILLE_CHARAC_ILLEGAL] = {'+', '-', '*', '[', ']', '(', ')'};
   premiere = li.find_first_of("\"") + 1;
   derniere = li.find_last_of("\"");
   string test(li.substr(premiere, derniere - premiere));
   // il faut rechercher et eliminer tout les caratÃšres spÃ©ciaux qu'ensight n'aime pas dans le nom
   for (j = 1; j < TAILLE_CHARAC_ILLEGAL; j++)
   {
      for (i = 0; i < test.size(); i++)
      {
         premiere = test.find_first_of(illegal[j]);
         if (premiere != -1) test.erase(premiere, 1);
      }
   }
   return test;
}

bool xExportTranslate::isSeparateur(char a) { return a == ','; }

bool xExportTranslate::isFinAc(char a) { return a == ')'; }

// ****** Partie traitement variables FICHIER ASCII

void xExportTranslate::traiteVariableScalaire(const string& ligne)
{
   int b = ligne.find_first_of("}");
   int j = 1;
   float hj;
   std::istringstream ins;
   string x("");
   for (; j <= b; j++)
   {
      if (!isSeparateur(ligne[j]) && !(ligne[j] == '}'))
      {
         x += ligne[j];
      }
      else
      {
         // if((x.size()>11)){

         ins.str(x);
         ins >> hj;
         fig.addScalar(hj);
         ins.clear();
         x.erase();
         //}
      }
   }
}

void xExportTranslate::traiteVariableVector(const string& ligne)
{
   int b = ligne.find_first_of("}");
   int j = 1, cpt = -1;
   float hj[3];
   std::istringstream ins;
   string x("");
   for (; j <= b; j++)
   {
      if (!isSeparateur(ligne[j]) && !(ligne[j] == '}'))
      {
         x += ligne[j];
      }
      else
      {
         cpt++;
         ins.str(x);
         ins >> hj[cpt];

         if (cpt == 2)
         {
            fig.addVect(hj[0], hj[1], hj[2]);

            cpt = -1;
         }

         ins.clear();
         x.erase();
      }
   }
}
void xExportTranslate::traiteVariableTensor(const string& ligne)
{
   int b = ligne.find_first_of("}");
   int j = 1, cpt = -1;
   float hj[9];
   float bte;
   std::istringstream ins;
   string x("");
   for (; j <= b; j++)
   {
      if (!isSeparateur(ligne[j]) && !(ligne[j] == '}'))
      {
         x += ligne[j];
      }
      else
      {
         cpt++;
         ins.str(x);
         ins >> bte;
         hj[cpt] = bte;
         // cout<<x<<"----------------------"<<hj[cpt]<<"       ********** "<<bte<<endl;
         if (cpt == 8)
         {
            fig.addTensor(hj[0], hj[1], hj[2], hj[3], hj[4], hj[5], hj[6], hj[7], hj[8]);
            // listeTens::const_iterator i;
            // i = fig.tensor.begin();
            // i->afficheTensor();
            // cout<<endl<<i->getTensor(0)<<endl;
            // exit(1);
            cpt = -1;
         }

         ins.clear();
         x.erase();
      }
   }
}

void xExportTranslate::traiteVariable(string ligne)
{
   // faire la distinction avec Scalar, Vector et Tensor.
   switch (fig.getTypeVariable())
   {
      case 1:
         traiteVariableScalaire(ligne);
         break;
      case 2:
         traiteVariableVector(ligne);
         break;
      case 3:
         traiteVariableTensor(ligne);
         break;
   }
}

void xExportTranslate::readGmshBinary()
{
   string b;
   int bufInt;
   int tab[29];
   int j = 0;
   double bufDbl;
   std::istringstream ins;
   getline(fpos_, b);
   // cout<<b<<endl;
   getline(fpos_, b);
   // cout<<b<<endl;
   getline(fpos_, b);
   // cout<<b<<endl;
   getline(fpos_, b);
   // la on vien de faire $view

   getline(fpos_, b);
   //  traitement de la ligne qui nous renseignement sur la suite que l'on a.
   bufInt = b.find_first_of(' ');
   b.substr(0, bufInt);
   //  std::cout<<"nom element : "<<b<<std::endl;
   // fig.addNomFigure(b.substr(0, bufInt));
   fig.addNomFigure(project_name_);
   string s = b.substr(bufInt + 1);
   // faire un swt=itch de 29 truc....
   int i;
   ins.str(s);

   for (i = 0; i < 29; i++)
   {
      ins >> tab[i];
      // std::cout<<tab[i];
   }
   for (i = 1; i < 29; i++)
   {
      if (tab[i] != 0)
      {
         readPosBinary(tab[i], i);  // on enregistre direct les element.
         j = i;
         // std::cout<<std::endl<<"nombre : "<<tab[i]<<" indice : "<<j<<std::endl;
      }
   }

   // ici on est renseigné sur le type a rentré donc c'est bon
   // on doit lire le pas de temps et la liste des valeurs, toujours à 1 donc :
   fpos_.read((char*)&bufInt, sizeof(int));
   fpos_.read((char*)&bufDbl, sizeof(double));
   // on doit maintenant lire la liste de double.
   // calcul du nombre de float
   int nbrNoeud = tab[j];

   i = fig.getTypeElem();
   if (i == 5) i--;
   bufInt = i * 3;  // nbre de coordonnée

   i = fig.getTypeVariable();
   if (i == 2) bufInt *= 2;
   if (i == 3) bufInt *= 4;

   double* lDbl = new double[nbrNoeud * bufInt];
   fpos_.read((char*)lDbl, (sizeof(double)) * nbrNoeud * bufInt);
   // fonction pour bien initialiser l'objet...
   for (i = 0; i < nbrNoeud; i++)
   {
      traiteElemBinaire(lDbl, i * bufInt);
   }
   delete[] lDbl;
   // faut voir comment recommencer la dedans.
}

void xExportTranslate::traiteElemBinaire(double* ptr, int cpt)
{
   int x = 0, y = 1, z = 2;
   int typeElem = fig.getTypeElem();
   int typeVar = fig.getTypeVariable();
   if (typeElem == 5) typeElem--;
   y *= typeElem;
   z *= typeElem;
   int i;
   for (i = 0; i < typeElem; i++)
   {
      fig.addElem(ptr[i + cpt + x], ptr[i + cpt + y], ptr[i + cpt + z]);

      if (typeVar == 1) fig.addScalar(ptr[i + (typeElem * 3)]);
      if (typeVar == 2)
         fig.addVect(ptr[(i * 3) + (typeElem * 3)], ptr[(i * 3) + (typeElem * 3) + 1], ptr[(i * 3) + (typeElem * 3) + 2]);
      if (typeVar == 3)
         fig.addTensor(ptr[(i * 3) + (typeElem * 3)], ptr[(i * 3) + (typeElem * 3) + 1], ptr[(i * 3) + (typeElem * 3) + 2],
                       ptr[(i * 3) + (typeElem * 3) + 3], ptr[(i * 3) + (typeElem * 3) + 4], ptr[(i * 3) + (typeElem * 3) + 5],
                       ptr[(i * 3) + (typeElem * 3) + 6], ptr[(i * 3) + (typeElem * 3) + 7], ptr[(i * 3) + (typeElem * 3) + 8]);
   }
}

void xExportTranslate::readPosBinary(int nbr, int cpt)
{
   switch (cpt)
   {
      case 1:
         firstOn(1, 1);
         break;
      case 2:
         firstOn(2, 1);
         break;
      case 3:
         firstOn(3, 1);
         break;
      case 4:
         firstOn(1, 2);
         break;
      case 5:
         firstOn(2, 2);
         break;
      case 6:
         firstOn(3, 2);
         break;
      case 7:
         firstOn(1, 3);
         break;
      case 8:
         firstOn(2, 3);
         break;
      case 9:
         firstOn(3, 3);
         break;
      case 10:
         firstOn(1, 5);
         break;
      case 11:
         firstOn(2, 5);
         break;
      case 12:
         firstOn(3, 5);
         break;
      case 13:
         firstOn(1, 4);
         break;
      case 14:
         firstOn(2, 4);
         break;
      case 15:
         firstOn(3, 4);
         break;
      case 16:
         firstOn(1, 8);
         break;
      case 17:
         firstOn(2, 8);
         break;
      case 18:
         firstOn(3, 8);
         break;
         // je n'ai pas implementé plus : pas de cube ou prism...
   }
}
}  // namespace xexport
