//
//  Listas.h
//  LLP
//
//  Created by cwilliam on 03/09/15.
//  Copyright (c) 2015 cwilliam. All rights reserved.
//

#ifndef __LLP__Listas__
#define __LLP__Listas__

#include <stdio.h>

typedef struct nodo_ {
    int dato;
    struct nodo_ *sig;
} nodo;


/*
    Funciones para Listas Ligadas.
 */

nodo *crea_nodo(int);
int  Enqueue(nodo **, nodo **, int);
int Dequeue(nodo **, nodo **, int *);
void Show_queue(nodo *);
void Realese_queue(nodo **);
int  Queue_size(nodo *);
void Initialization(nodo **);

#endif /* defined(__LLP__Listas__) */
