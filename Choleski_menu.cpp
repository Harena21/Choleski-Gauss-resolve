/*
 * Résolution de système linéaire par Méthode Gauss et la Méthode de Choleski
 * V 1.0
 * 2021-06-22
 * Auteur(s): HARENA Antenaina Lucka
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <locale>
#include <iomanip>
#include <string>
using namespace std;

/// General helpers
float **newMat(int rows, int cols);
template <class T>
T *newVect(size_t dim);
void displayMat(size_t dim, float **A);
void displayVec(size_t dim, float *v);

class Lsolver
{
public:
    Lsolver(string filename);
    ~Lsolver();
    void displayResultGauss();      //pour afficher spécifiquement le résultat de la méthode de Gauss
    void displayResultChole();      //pour afficher spécifiquement le résultat de la méthode de Choleski
    void gaussElim();               //pour résoudre l'équation par l'élimination de Gauss
    void cholesky();                //Traiter le système par la méthode de Choleski
    void solveTriangSupGauss();     //pour résoudre la matrice triangulaire supérieure eu dans la méthode de Gauss 
    void solveTriangSupChole();     //pour résoudre la matrice triangulaire supérieure eu dans la méthode de Choleski
    void solveTriangInf();          //pour résoudre la matrice triangulaire inférieure dans la méthode de Choleski 
    void transposeMat();
    size_t getdim() { return dim; }
    float **getMatA() { return A; }
    float **getMatB() { return B; }     //retourne la matrice B pour la nouvelle méthode
    float *getY(){ return y;}           //retourne le vecteur Y issu de la méthode de Choleski
    float *getX(){ return x;}
    float *getRhs() { return b; }

private:
    size_t dim;
    float **A; // matrice du problème A.x=b
    float **B; //nouvelle matrice pour traiter le système linéaire dans le cas de la méthode de Choleski
    float *b, *x, *y; // second membre et inconnu du problème
    int *p;           // pointeur utilisé pour les permutations de ligne
};
class Menu      //nouvelle classe pour afficher un menu dans le main
{
    public:
    Menu();     //constructeur de la classe Menu
    int display();      //méthode pour afficher le menu 

    private:
    string tab[3];      
    int choice;
};

Lsolver::Lsolver(string filename)
{
    /// ouverture du fichier de donn�es
    size_t i(0), j(0);
    ifstream fichier(filename, ios::in);
    if (fichier)
    {
        fichier >> dim;
        /// allocation de la matrice, du second membre et de la solution
        A = newMat(dim, dim);
        B = newMat(dim, dim);
        b = newVect<float>(dim);
        x = newVect<float>(dim);
        y = newVect<float>(dim);
        p = newVect<int>(dim); // pointeur pour la permutation de lignes
                               /// remplissage des donn�es
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
            {
                fichier >> A[i][j];
            }
        }
        for (i = 0; i < dim; i++)
        {
            fichier >> b[i];
        }
        for (i = 0; i < dim; i++)
        {
            p[i] = i;
            x[i] = 0;
        }
        fichier.close();
    }
    else
    {
        cout << "Donn�es non trouv�es..." << endl;
    }
}

Lsolver::~Lsolver()
{
    delete[] p;
    delete[] x;
    delete[] y;
    delete[] b;
    for (size_t i = 0; i < dim; i++)
        delete[] A[i];
    delete[] A;
    for (size_t i = 0; i < dim; i++)
        delete[] B[i];
    delete[] B;
}

void Lsolver::displayResultGauss()  //Afficher le résultat après la méthode de Gauss
{
    float eps(1e-6);
    cout << "\nLe probl�me triangularis�:" << endl;
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < dim - 1; j++)
        {
            if (fabs(A[p[i]][j]) > eps || j >= i)
                cout << A[p[i]][j] << ".x" << j + 1 << " + ";
            else
                cout << "       ";
        }
        cout << A[p[i]][dim - 1] << ".x" << dim << " = " << b[p[i]] << endl;
    }
    cout << "\nLa solution:" << endl;
    for (size_t i = 0; i < dim; i++)
        cout << "x" << i + 1 << " = " << x[i] << endl;
}

void Lsolver::gaussElim(){
    size_t t(0), lpiv(0);    /// ligne de pivot courant
    float piv(0);
    for(size_t k=0; k<dim; k++){

/// Recherche du plus grand pivot
        lpiv = k;
        piv = fabs(A[p[k]][k]);
        for(size_t i=k+1; i<dim; i++){
            if(fabs(A[p[i]][k]) > piv){
                lpiv = i;
                piv = fabs(A[p[i]][k]);
            }
        }

/// Permutation de lignes
        t = p[k];
        p[k] = p[lpiv];
        p[lpiv] = t;

/// Elimination sous la ligne de pivot
        for(size_t i=k+1; i<dim; i++){
			A[p[i]][k] /= A[p[k]][k];
            for(size_t j=k+1; j<dim; j++)
                A[p[i]][j] -= (A[p[i]][k] * A[p[k]][j]);
            b[p[i]] -= (A[p[i]][k] * b[p[k]]);
            A[p[i]][k] = 0;
        }
    }
}
void Lsolver::solveTriangSupGauss()
{
    float s(0);
    int i(0), j(0);
    for (i = dim - 1; i >= 0; i--)
    { /// Must go backward
        for (j = i + 1, s = 0; j < int(dim); j++)
            s += (A[p[i]][j] * x[j]);
        x[i] = (b[p[i]] - s) / A[p[i]][i];
    }
}
void Lsolver::cholesky()    //traitement du système par la méthode de Choleski
{
    size_t s(0);
    B[0][0] = sqrt(A[0][0]);        //Toujours vrai donc initialisé à l'extérieur de la boucle

    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j <= i; j++)
        {
            s = 0;
            if (j == 0)         
            {
                B[i][j] = A[i][j] / B[0][0];
            }
            else
            {
                for (size_t k = 0; k < j; k++)
                {
                    s += B[i][k] * B[j][k];
                }
                if (i == j)
                {
                    B[j][j] = sqrt(A[j][j] - s);
                }
                else
                {
                    B[i][j] = (A[i][j] - s) / B[j][j];
                }
            }
        }
    }
}

void Lsolver::solveTriangSupChole()
{
    float s(0);
    for(int i=dim-1 ; i>=0 ; i--)       //pour trouver les solutions x de l'équation A.X=B par la méthode de Choleski
    {
        s = 0;
        for(int j=i+1 ; j<dim ; j++)
        {
            s += B[i][j]*x[j];
        }
        x[i] = (y[i] - s)*(1/B[i][i]);
    }
}
void Lsolver::solveTriangInf()
{
    float s(0);
    y[0] = b[0] / B[0][0];
    for (int i = 1; i < dim; i++)   //pour trouver les solutions y de l'équation B.Y=b vu dans la méthode de Cholesksi
    {
        s = 0;
        for (int j = 0; j < i; j++)
        {
            s += (B[i][j] * y[j]);
        }
        y[i] = (b[i] - s) / B[i][i];
    }
}
void Lsolver::transposeMat()
{
    float tmp(0);
    for (int i = 0; i < dim; i++) // pour transposer la matrice triangle inférieure en matrice triangle supérieure
    {
        for (int j = i+1; j < dim; j++)
        {
            tmp = B[i][j];
            B[i][j] = B[j][i];
            B[j][i] = tmp;
        }
    }
}

Menu::Menu()
{
    tab[0] = "Résoudre par la méthode de Gauss";        //Tous les choix possibles du menu
    tab[1] = "Résoudre par la méthode de Choleski";
    tab[2]= "Quitter";
    choice = 0;
}
int Menu::display()
{
    for(int i=0; i<2; i++)
    {
        cout << i+1 << "\t" << tab[i] << endl;
    }
    cout << "0" << "\t" << tab[2] << endl;
    cout <<"Entrez votre choix : " ; cin >> choice; 
    return choice;                              // retourne le choix de l'utilisateur
}
void Lsolver::displayResultChole()       //pour afficher le résultat de la méthode de Choleski dans main
{
    cout << "Le problème triangularisé : " << endl;
    for(int i=0 ; i<dim ; i++)
    {
        for(int j=i ; j<dim ; j++)
        {
            cout << B[i][j] << "." << "x" << j+1 ;
            if(j<dim-1)
            {
                cout << " + ";
            }
            else
            {
                cout << " = ";
            }
        }
         cout << y[i] <<endl;
    }

    cout << "\nLa solution : " << endl;
    for(int i=0 ; i<dim ; i++)
    {
        cout << "x" << i+1 << " = " << x[i] << endl;
    }
}

int main()
{
    cout << "Locale: " << setlocale(LC_ALL, "") << endl;
    cout << "Résolution d'un système linéaire A.x=b" << endl;

    /// Données
    Lsolver solver("texte");
    Menu menu;
    int choice(0);

    cout << "Voici la matrice A :" << endl ;
    displayMat(solver.getdim(),solver.getMatA());
    cout << "Le second membre :" << endl;
    displayVec(solver.getdim(),solver.getRhs()); 

    /// Traitement et Affichage 
        choice = menu.display();

        if(choice == 1)
        {
            solver.gaussElim();
            solver.solveTriangSupGauss();
            solver.displayResultGauss();
        }
        else if(choice == 2)
        {
            cout << "On suppose que la matrice associée au sytème est symétrique et définie positive :" << endl;
            solver.cholesky();
            cout << "La matrice après la méthode de Choleski:" << endl;
            displayMat(solver.getdim(),solver.getMatB());
            solver.solveTriangInf();

            solver.transposeMat();
            cout << "La Matrice transposée : " << endl;
            displayMat(solver.getdim(),solver.getMatB());
            cout<< "Le second membre qui devient :" << endl;
            displayVec(solver.getdim(),solver.getY());
            
            solver.solveTriangSupChole();

            solver.displayResultChole();
        }
      
    return 0;
}

void displayMat(size_t dim, float **A)
{
    cout << "[";
    for (size_t i = 0; i < dim; i++)
    {
        if (i == 0)
            cout << "[";
        else
            cout << " [";
        for (size_t j = 0; j < dim; j++)
        {
            if (j < dim - 1)
                cout << setw(8) << setprecision(5) << A[i][j] << "  ";
            else
                cout << setw(8) << setprecision(5) << A[i][j];
        }
        if (i < dim - 1)
            cout << "]" << endl;
        else
            cout << "]]" << endl;
    }
}
void displayVec(size_t dim, float *v){
	for(size_t i=0; i<dim; i++){
		cout << " [" << v[i] << "]" << endl;
	}
}
template <class T>
T *newVect(size_t dim)
{
    T *v(NULL);
    v = new T[dim]; // Qu'est-ce que je fais en cas d'erreur
    return v;
}

float **newMat(int rows, int cols)
{
    float **mat(NULL);
    mat = newVect<float *>(rows);

    for (int i = 0; i < rows; i++)
        mat[i] = newVect<float>(cols);
    return mat;
}

