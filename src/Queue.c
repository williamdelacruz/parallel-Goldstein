//
//  Queue
//
//  Created by cwilliam on 03/09/15.
//  Copyright (c) 2015 cwilliam. All rights reserved.
//

#include "Queue.h"
#include <stdlib.h>

/*
    Funcion para crear un nodo
 */
nodo *crea_nodo(int dato) {
    nodo *nuevo;
    
    nuevo = (nodo *)malloc(sizeof(nodo));
    if (nuevo==NULL) {
        printf("No hay memoria.\n");
        return NULL;
    }
    nuevo->dato=dato;
    nuevo->sig=NULL;
    return nuevo;
}


/* 
    Inserta al inicio
 */
/*
int Enqueue(nodo **lista, int info) {
    nodo *nuevo;
    
    nuevo = crea_nodo(info);
    if (nuevo==NULL)
        return 0;
    
    if (*lista==NULL) {
        *lista = nuevo;
    }
    else {
        nuevo->sig = *lista;
        *lista = nuevo;
    }
    return 1;
}
*/


int Enqueue(nodo **Q, nodo **endQ, int info)
{
    nodo *nuevo;
    nuevo = crea_nodo(info);

    if (*Q==NULL) {
        *Q = nuevo;
        *endQ = nuevo;
    }
    else {
    	(*endQ)->sig=nuevo;
    	*endQ = nuevo;
    }
    return 1;
}


/*
    Elinina al final
 */
int Dequeue(nodo **Q, nodo **endQ, int *dato)
{
	nodo *tmp;
    if (*Q==NULL)
        return -1;
    else {
    	*dato = (*Q)->dato;
        tmp = *Q;
        *Q = (*Q)->sig;
        free(tmp);
        
        if (*Q==NULL)
        	*endQ = NULL;
    }
    
    return 0;
}


/*
    Recorre e imprime una lista ligada.
 */
void Show_queue(nodo *lista) {
    nodo *ptr;
    
    printf("\n");
    for (ptr=lista; ptr!=NULL; ptr=ptr->sig)
        printf("%d\t",ptr->dato);
    printf("\n");
}


/*
    Borra los nodos de una lista ligada.
 */
void Realese_queue(nodo **lista) {
    nodo *ptr, *tmp;
    
    ptr = *lista;
    while (ptr!=NULL) {
        tmp = ptr;
        ptr = ptr->sig;
        free(tmp);
    }
    *lista = NULL;
}

/*
    Cuenta el numero de nodos de una lista ligada.
 */
int Queue_size(nodo *lista) {
    nodo *ptr;
    int k=0;
    
    for (ptr=lista; ptr!=NULL; ptr=ptr->sig, k++);
    return k;
}



/*
    Inicializa a NULL un apuntador a una lista ligada.
 */
void Initialization(nodo **lista) {

    *lista = NULL;
}

