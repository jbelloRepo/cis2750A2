#include "mol.h"
#define PI 3.14159265358979323846 // define PI
#define DEBUG_OFF				  // for debugging

/************************************
******** Function Definition ********
************************************/

int compare_atom_ptr(const void *a, const void *b);
int compare_bond_ptr(const void *a, const void *b);

/******************************************
		Functions for the atom
******************************************/

void atomset(atom *atom, char element[3], double *x, double *y, double *z) // setter for atom
{
	//* Copy elements into struct members
	strcpy(atom->element, element); // strcpy to struct
	atom->x = *x;
	atom->y = *y;
	atom->z = *z;

#ifdef DEBUG_
	printf("The string read from char array is: %s \n", atom->element);
#endif
}

void atomget(atom *atom, char element[3], double *x, double *y, double *z)
{
	//* Copy elements into struct members
	if (atom != NULL)
	{
		strcpy(element, atom->element); // strcpy to struct
		*x = atom->x;
		*y = atom->y;
		*z = atom->z;
	}
}

/******************************************
		Functions for the molecule
******************************************/

molecule *molmalloc(unsigned short atom_max, unsigned short bond_max) // dym. allocate space for the molecule -> return pointer of type malloc
{
	molecule *molecule_ptr = (molecule *)malloc(sizeof(molecule)); //* allocate molecule on the heap

#ifdef DEBUG_
	printf("address of molecule_ptr in header file : %p\n", (void *)molecule_ptr);
#endif

	if (molecule_ptr != NULL) //* allocate space for "atoms" array & atom_ptrs
	{
		molecule_ptr->atom_max = atom_max;
		molecule_ptr->atom_no = 0;									   // set to 0 initially
		molecule_ptr->atoms = (atom *)malloc(atom_max * sizeof(atom)); // unsigned short vs int

#ifdef DEBUG_
		printf("address of molecule_ptr->atoms in header file : %p\n", (void *)molecule_ptr->atoms); //! error checking
#endif
		molecule_ptr->atom_ptrs = (atom **)malloc(atom_max * sizeof(atom *));
		molecule_ptr->bond_max = bond_max;
		molecule_ptr->bonds = (bond *)malloc(bond_max * sizeof(bond));
		molecule_ptr->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
		molecule_ptr->bond_no = 0; // set to 0 initially
	}
	else
	{
		return NULL;
	}
	// molecule *mol_malloced_address = (molecule *)malloc(sizeof(molecule)); // pointer to molecule
	return molecule_ptr;
}

//! Add NULL for malloc if there is no space FIXME: Add NUL
void molappend_atom(molecule *molecule, atom *atom) // FIXME: after re-allocating you have to reassign the atom pointer, same for the bonds
{
	struct atom *a1, **a2; //! TEMP VARS

	if (molecule->atom_no == molecule->atom_max) //* check atom_no & atom_max before appending to array
	{
		if (molecule->atom_max == 0) //! TEST THIS when atom_max = 0
		{
			molecule->atom_max += 1; //* add 1 if atom_max = 0
		}
		else
		{
			molecule->atom_max *= 2; //* double size of atom_max
		}

		/* Reallocating Memory (ATOM) */

#ifdef DEBUG_ON
		printf("\n==================================== [MOL.H] This is for molappend_atom() ==================================\n");
		printf("The (ATOM) array is now full. REALLOCATING!!!....\n"); //! error checking
#endif
																	   // *(molecule->atoms) = (atom *)realloc(molecule->atoms, molecule->atom_max);
		a1 = (struct atom *)realloc(molecule->atoms, molecule->atom_max * sizeof(struct atom)); // you cannot use atom* here

#ifdef DEBUG_ON
		printf("Address pointed to by {molecule->atoms} (old memory address): %p\n", (void *)molecule->atoms);	  //! error checking
		printf("Address pointed to by {temp a1} in header file (new memory address #REALLOC): %p\n", (void *)a1); //! error checking
#endif

		if (a1 == NULL) //! add more statements or return 0 --CHECK THIS
		{
			exit(EXIT_FAILURE);
#ifdef DEBUG_ON
			printf("No additional heap for atoms *atom \n");
#endif
		}

		else
		{

			molecule->atoms = a1;

#ifdef DEBUG_ON
			printf("Address assigned by {temp a1} to {molecule->atoms} in header file (new memory address #REALLOC): %p\n", (void *)molecule->atoms); //! error checking
#endif
		}

		a2 = (struct atom **)realloc(molecule->atom_ptrs, molecule->atom_max * sizeof(struct atom *));

#ifdef DEBUG_ON
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("Address pointed to by {molecule->atom_ptr} (old memory address): %p\n", (void *)molecule->atom_ptrs); //! error checking
		printf("Address pointed to by {temp a2} in header file (new memory address #REALLOC) : %p\n", (void *)a2);	  //! error checking
#endif
		if (a2 == NULL) //! add more statements or return 0 --CHECK THIS
		{
			exit(EXIT_FAILURE);

#ifdef DEBUG_ON
			printf("No additional heap for atom **atoms \n");
#endif
		}
		else
		{
			molecule->atom_ptrs = a2;

#ifdef DEBUG_ON
			printf("Address assigned by {temp a2} to {molecule->atom_ptrs} in header file(new memory address #REALLOC) : %p\n", (void *)molecule->atom_ptrs); //! error checking
#endif
		}

		for (int i = 0; i < molecule->atom_no; i++)
		{
			molecule->atom_ptrs[i] = &molecule->atoms[i];
		}

#ifdef DEBUG_ON
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("Address pointed of index[0] in atoms (after re-alloc): %p\n", (void *)&molecule->atoms[0]);			 //! error checking
		printf("Address pointed to by {molecule->atom_ptr} (after re-alloc): %p\n", (void *)molecule->atom_ptrs[0]); //! error checking
																													 // printf("Address pointed of index[1] in atoms (after re-alloc): %p\n", (void *)&molecule->atoms[1]);			 //! error checking
																													 // printf("Address pointed to by {molecule->atom_ptr} (after re-alloc): %p\n", (void *)molecule->atom_ptrs[1]); //! error checking
#endif

		// free(a1); //? FREE THIS ?
		// free(a2);
	}

	if ((molecule->atom_no > 0) && (molecule->atom_no < molecule->atom_max)) //* check atom_no is greater than 0 and lesser than atom_max
	{
		strcpy(molecule->atoms[molecule->atom_no].element, atom->element);
		molecule->atoms[molecule->atom_no].x = atom->x;
		molecule->atoms[molecule->atom_no].y = atom->y;
		molecule->atoms[molecule->atom_no].z = atom->z;
		molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[molecule->atom_no]; // pointer to the first atom
		molecule->atom_no += 1;														  // increment after addition
	}

	if ((molecule->atom_no == 0) && (molecule->atom_max != 0)) //* check if the molecule is initially empty ( executes for first element only)
	{
		//! index is static for first element in array[0]
		strcpy(molecule->atoms[0].element, atom->element);
		molecule->atoms[0].x = atom->x;
		molecule->atoms[0].y = atom->y;
		molecule->atoms[0].z = atom->z;
		molecule->atom_ptrs[0] = &molecule->atoms[0]; // pointer to the first atom
		molecule->atom_no += 1;						  // increment after adding to atoms

#ifdef DEBUG_ON
		printf("address pointed to by molecule->atom_ptrs[0] in header file : %p\n", (void *)molecule->atom_ptrs[0]); //! error checking
		printf("The value of atom_max is %d \n", molecule->atom_max);
		printf("The number of bonds currently in the molecule is: %d \n", molecule->bond_no);
		printf("===============================================================================================================\n");
#endif
	}
}

/******************************************
		Functions for the bond
*******************************************/

void bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs) //* setter for the bond variable
{
	bond->a1 = a1;
	bond->a2 = a2;

#ifdef DEBUG_ON
	printf("\n==================================== [MOL.C] This is for bondset() =======================================\n");
	printf("Address pointed to by bond->a1 in the header file: %p\n", (void *)bond->a1); //! error checking
	printf("Address pointed to by bond->a2 in the header file: %p\n", (void *)bond->a2); //! error checking
#endif

	bond->epairs = epairs;
}

void bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs)
{
	*a1 = (bond->a1); //*a1 yields pointer to atom; (bond->a1) returns pointer to atom
	*a2 = (bond->a2);
	*epairs = bond->epairs;
}

void molappend_bond(molecule *molecule, bond *bond) // appends the bonds
{													//! bond_no initially set to 0

	struct bond *b1, **b2; // TEMP variables

#ifdef DEBUG_
	printf("\n==================================== [MOL.H] This is for molappend_bond() =======================================\n");
	printf("The value of bond_no and bond_max is: %d , %d\n", molecule->bond_no, molecule->bond_max); //! error checking
#endif

	if (molecule->bond_no == molecule->bond_max) //* check if bond_no and bond_max are equal
	{
		if (molecule->bond_max == 0) //! TEST THIS when bond_max = 0
		{
			molecule->bond_max += 1; //* add 1 if bond_max = 0
		}
		else
		{
			molecule->bond_max *= 2; //* double size of bond_max
		}

		/* Reallocating Memory (BOND) */

#ifdef DEBUG_
		printf("\n==================================== [MOL.H] This is for molappend_bond() =======================================\n");
		printf("The (BOND) array is now full. REALLOCATING!!!....\n"); //! error checking
#endif

		b1 = (struct bond *)realloc(molecule->bonds, molecule->bond_max * sizeof(struct bond)); // you cannot use atom* here

#ifdef DEBUG_ON
		printf("Address pointed to by {molecule->bonds} (old memory address): %p\n", (void *)molecule->bonds);	  //! error checking
		printf("Address pointed to by {temp b1} in header file (new memory address #REALLOC): %p\n", (void *)b1); //! error checking
#endif

		if (b1 == NULL) //! add more statements or return 0 --CHECK THIS
		{
			exit(EXIT_FAILURE);

#ifdef DEBUG_ON
			printf("No additional heap for atoms *atom \n");
#endif
		}

		else
		{
			molecule->bonds = b1;

#ifdef DEBUG_ON
			printf("Address assigned by {temp b1} to {molecule->bonds} in header file (new memory address #REALLOC): %p\n", (void *)molecule->bonds); //! error checking
#endif
		}

		b2 = (struct bond **)realloc(molecule->bond_ptrs, molecule->bond_max * sizeof(struct bond *));

#ifdef DEBUG_ON
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("Address pointed to by {molecule->bond_ptr} (old memory address): %p\n", (void *)molecule->bond_ptrs); //! error checking
		printf("Address pointed to by {temp b2} in header file (new memory address #REALLOC): %p\n", (void *)b2);	  //! error checking
#endif

		if (b2 == NULL) //! add more statements or return 0 --CHECK THIS
		{
			exit(EXIT_FAILURE);

#ifdef DEBUG_ON
			printf("No additional heap for atom **atoms \n");
#endif
		}
		else
		{
			molecule->bond_ptrs = b2;

#ifdef DEBUG_ON
			printf("Address assigned by {temp b2} to molecule->bond_ptrs in header file(new memory address #REALLOC) : %p \n\n", (void *)molecule->bond_ptrs); //! error checking
#endif
		}

		for (int i = 0; i < molecule->bond_no; i++)
		{
			molecule->bond_ptrs[i] = &molecule->bonds[i];
		}

#ifdef DEBUG_ON
		printf("-----------------------------------------------------------------------------------------------------------\n");
		printf("Address pointed of index[0] in atoms (after re-alloc): %p\n", (void *)&molecule->bonds[0]);			 //! error checking
		printf("Address pointed to by {molecule->atom_ptr} (after re-alloc): %p\n", (void *)molecule->bond_ptrs[0]); //! error checking
#endif
	}

	if ((molecule->bond_no > 0) && (molecule->bond_no < molecule->bond_max)) //* check bond_no is greater than 0 and lesser than bond_max
	{
#ifdef DEBUG_
		printf("\n==================================== [MOL.H] This is for molappend_bond() =======================================\n");
		printf("EXE: When bond_no is greater than 0 and less than BOND MAX \n"); //! error checking
		printf("The value of bond_no and bond_max is: %d , %d\n", molecule->bond_no, molecule->bond_max);
#endif
		molecule->bonds[molecule->bond_no].a1 = bond->a1; //* assign the first member from parameter passed to molecule
		molecule->bonds[molecule->bond_no].a2 = bond->a2;
		molecule->bonds[molecule->bond_no].epairs = bond->epairs;
		molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[molecule->bond_no]; //! CHECK THIS
		molecule->bond_no += 1;														  // increment after addition
	}

	if ((molecule->bond_no == 0) && (molecule->bond_max != 0)) //* check if the molecule is initially empty (executes for first element only)
	{
//! index is static for first element in array[0]
#ifdef DEBUG_
		printf("\n==================================== [MOL.H] This is for molappend_bond() =======================================\n");
		printf("EXE: When there is no bond in the BOND ARRAY --> Append bond to index [0] in the array \n"); //! error checking
#endif
		molecule->bonds[0].a1 = bond->a1; //* assign the first member from parameter passed to molecule
		molecule->bonds[0].a2 = bond->a2;
		molecule->bonds[0].epairs = bond->epairs;
		molecule->bond_ptrs[0] = &molecule->bonds[0]; //! CHECK THIS
		molecule->bond_no += 1;

#ifdef DEBUG_ON
		printf("Address pointed to by molecule->bond_ptrs[0] in header file : %p\n", (void *)molecule->bond_ptrs[0]); //! error checking
		printf("The value of bond_max is %d \n", molecule->bond_max);
		printf("The number of bonds currently in the molecule is: %d \n", molecule->bond_no);
		// printf("===========================================================================================================\n");
#endif
	}
}

/******************************************
		Additional functions
*******************************************/

int compare_atom_ptr(const void *a, const void *b) // compare atom pointers
{
	// pass in an array of type pointer to atom

	double double_a, double_b;
	atom *atom_a, *atom_b;

	atom_a = *(atom **)(a); // de-reference a to get the atom pointed to by the pointer
	double_a = atom_a->z;	// get the value pointed to by a (which is the z value)

	atom_b = *(atom **)(b);
	double_b = atom_b->z;

#ifdef DEBUG_
	printf("The value of a is: %lf \n", double_a);
	printf("The value of b is: %lf \n", double_b);
#endif

	return ((int)(double_a - double_b));
}

int compare_bond_ptr(const void *a, const void *b)
{
	double double_a, double_b;
	bond *bond_a, *bond_b;

	bond_a = *(bond **)(a);							// de-reference a to get the bond pointed to by the pointer
	double_a = (bond_a->a1->z + bond_a->a2->z) / 2; // get the value pointed to by a (which is the z value)

	bond_b = *(bond **)(b);
	double_b = (bond_b->a1->z + bond_b->a2->z) / 2;

#ifdef DEBUG_
	printf("The value of a is: %lf \n", double_a);
	printf("The value of b is: %lf \n", double_b);
#endif

	return ((int)(double_a - double_b));
}

void molsort(molecule *molecule) // This function sorts an array of pointers (pointing to atoms in an array)
{
#ifdef DEBUG_ON
	printf("\n==================================== [MOL.C] This is for molsort() [atom_ptrs] ==============================\n");
	printf("The number of atoms in the molecule is %d\n", molecule->atom_no);
	printf("BEFORE atom_ptr is sorted:");
	for (int i = 0; i < molecule->atom_no; i++)
	{
		if (i == 0)
		{
			printf(" %p", (void *)molecule->atom_ptrs[i]);
		}
		else
		{
			printf(" %p", (void *)molecule->atom_ptrs[i]);
		}
	}
	printf("\n");
#endif

	qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), compare_atom_ptr); // sort the atom_ptrs

#ifdef DEBUG_ON
	printf("AFTER  atom_ptr is sorted:");
	for (int i = 0; i < molecule->atom_no; i++)
	{
		if (i == 0)
		{
			printf(" %p", (void *)molecule->atom_ptrs[i]);
		}
		else
		{
			printf(" %p", (void *)molecule->atom_ptrs[i]);
		}
	}
	printf("\n");
#endif

#ifdef DEBUG_ON
	printf("\n==================================== [MOL.C] This is for molsort() [bond_ptrs] ==============================\n");
	printf("The number of bonds in the molecule is %d\n", molecule->bond_no);
	printf("BEFORE bond_ptr is sorted:");
	for (int i = 0; i < molecule->bond_no; i++)
	{
		if (i == 0)
		{
			printf(" %p", (void *)molecule->bond_ptrs[i]);
		}
		else
		{
			printf(" %p", (void *)molecule->bond_ptrs[i]);
		}
	}
	printf("\n");
#endif

	qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), compare_bond_ptr); // sort the bond_ptrs

#ifdef DEBUG_ON
	printf("AFTER  bond_ptr is sorted:");
	for (int i = 0; i < molecule->bond_no; i++)
	{
		if (i == 0)
		{
			printf(" %p", (void *)molecule->bond_ptrs[i]);
		}
		else
		{
			printf(" %p", (void *)molecule->bond_ptrs[i]);
		}
	}
	printf("\n");
#endif
}

void xrotation(xform_matrix xform_matrix, unsigned short deg)
{
	//* the xform_matrix act as a pointer since an array is passed

	//  ---           ---
	// | 1,    0,    0   |
	// | 0, cos(), -sin()|
	// | 0, sin(), cos() |
	//  ---           ---

	double rad = deg * (PI / 180.0);
	xform_matrix[0][0] = 1;
	xform_matrix[0][1] = 0;
	xform_matrix[0][2] = 0;
	xform_matrix[1][0] = 0;
	xform_matrix[1][1] = cos(rad);
	xform_matrix[1][2] = -sin(rad);
	xform_matrix[2][0] = 0;
	xform_matrix[2][1] = sin(rad);
	xform_matrix[2][2] = cos(rad);
}
void yrotation(xform_matrix xform_matrix, unsigned short deg)
{
	//  ---           ---
	// | cos(), 0, sin() |
	// |    0   1,   0   |
	// |-sin(), 0, cos() |
	//  ---           ---

	double rad = deg * (PI / 180.0);
	xform_matrix[0][0] = cos(rad);
	xform_matrix[0][1] = 0;
	xform_matrix[0][2] = sin(rad);
	xform_matrix[1][0] = 0;
	xform_matrix[1][1] = 1;
	xform_matrix[1][2] = 0;
	xform_matrix[2][0] = -sin(rad);
	xform_matrix[2][1] = 0;
	xform_matrix[2][2] = cos(rad);
}
void zrotation(xform_matrix xform_matrix, unsigned short deg)
{
	//  ---             ---
	// | cos(), -sin(),  0 |
	// | sin()   cos(),  0 |
	// |   0,      0,    1 |
	//  ---           ---

	double rad = deg * (PI / 180.0);
	xform_matrix[0][0] = cos(rad);
	xform_matrix[0][1] = -sin(rad);
	xform_matrix[0][2] = 0;
	xform_matrix[1][0] = sin(rad);
	xform_matrix[1][1] = cos(rad);
	xform_matrix[1][2] = 0;
	xform_matrix[2][0] = 0;
	xform_matrix[2][1] = 0;
	xform_matrix[2][2] = 1;
}

void mol_xform(molecule *molecule, xform_matrix matrix)
{
	// apply matrix vector multiplication
	double x, y, z; // position of atom relative to molecule

#ifdef DEBUG_ON
	printf("\n==================================== [MOL.C] This is for mol_xform() ==============================\n");
	for (int i = 0; i < molecule->atom_no; i++)
	{
		printf("BFR MATRIX TRANSFORM: %lf, %lf, %lf\n", molecule->atoms[i].x, molecule->atoms[i].y, molecule->atoms[i].z);
	}
#endif

	for (int i = 0; i < molecule->atom_no; i++) // apply transformation to x,y,z values in molecule
	{
		x = molecule->atoms[i].x;
		y = molecule->atoms[i].y;
		z = molecule->atoms[i].z;

		molecule->atoms[i].x = (matrix[0][0] * x) + (matrix[0][1] * y) + (matrix[0][2] * z);
		molecule->atoms[i].y = (matrix[1][0] * x) + (matrix[1][1] * y) + (matrix[1][2] * z);
		molecule->atoms[i].z = (matrix[2][0] * x) + (matrix[2][1] * y) + (matrix[2][2] * z);
	}

#ifdef DEBUG_ON
	printf("\n==================================== [MOL.C] This is for mol_xform() ==============================\n");
	for (int i = 0; i < molecule->atom_no; i++)
	{
		printf("ATR MATRIX TRANSFORM: %lf, %lf, %lf\n", molecule->atoms[i].x, molecule->atoms[i].y, molecule->atoms[i].z);
	}
#endif
}

void molfree(molecule *ptr)
{
	free(ptr->atoms);
	free(ptr->atom_ptrs);
	free(ptr->bonds);
	free(ptr->bond_ptrs);
	free(ptr);
}

molecule *molcopy(molecule *src)
{
	// src is a molecule passed

	// // atom struct members
	// double x, y, z;
	// char element[3];

	// // bond struct members
	// atom *a1, *a2;
	// unsigned char epairs;

	atom anAtom;		// an atom variable
	bond aBond;			// a bond variable
	molecule *moleCopy; // molecule variable

	unsigned short atom_max = src->atom_max, bond_max = src->bond_max;

	// printf("\nThe atom max and bond max: %d %d\n", atom_max, bond_max);
	moleCopy = molmalloc(atom_max, bond_max); // dym allocate space for mole copy

	// set strcut members in molecule
	moleCopy->atom_max = src->atom_max;
	// moleCopy->atom_no = src->atom_no; //! FIXME -- do not set this
	moleCopy->bond_max = src->bond_max;
	// moleCopy->bond_no = src->bond_no; //! FIXME -- do not set this

	// append atoms in atom arrays
	for (int i = 0; i < src->atom_no; i++)
	{
		//? Create the atom at array index
		strcpy(anAtom.element, src->atoms[i].element);
		anAtom.x = src->atoms[i].x;
		anAtom.y = src->atoms[i].y;
		anAtom.z = src->atoms[i].z;

#ifdef DEBUG_ON
		printf("\n==MOLCOPY== Adding atom to the molecule copy\n");
#endif

		molappend_atom(moleCopy, &anAtom); // append all atoms to molcopy

#ifdef DEBUG_ON
		printf("==MOLCOPY== A new atom has been added. TOTAL COUNT: %d \n", i + 1);
#endif
	}

	for (int i = 0; i < src->bond_no; i++)
	{
		//? Create the bond at bond index
		aBond.epairs = src->bonds[i].epairs;
		aBond.a1 = src->bonds[i].a1;
		aBond.a2 = src->bonds[i].a2;

		molappend_bond(moleCopy, &aBond); // append all bonds
	}

	return moleCopy;
}
