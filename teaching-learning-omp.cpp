#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <random>

#define MAX_SUBJECT 5
#define MAX_NUMBER_SOLUTION 9
#define MIN_NUMBER_SOLUTION 0
#define MAX_STUDENT 15
#define NUM_GENERATIONS 100
#define MAX_EVALUATIONS 640

// Semilla aleatoria utilizando el reloj del sistema
std::random_device rd;

// Generador de números aleatorios
std::mt19937 gen(rd());

// Distribución uniforme en el rango [0, MAX_NUMBER_SOLUTION)
std::uniform_int_distribution<> dis(0, MAX_NUMBER_SOLUTION);

int num_evaluations = 0;
int num_evaluations_found_solution = 0;

using namespace std;

//type student
struct Student {
    int solution[MAX_SUBJECT];
    float fx;
};

// Imprimir una solución con su índice, valores y aptitud
void print_student(Student *student, int index) {
    cout << "Student " << index << ": ";
    for (int i = 0; i < MAX_SUBJECT; ++i) {
        cout << student->solution[i] << " ";
    }
    cout << "| fx: " << student->fx;
}

Student* select_the_teacher(Student classroom[], int max_students){
    Student* teacher = &classroom[0];
    for (int student = 0; student < max_students; student++)
        if(classroom[student].fx > teacher->fx) teacher = &classroom[student];
    return teacher; // the first student is the best in the classroom
}

/**
 * Funcion que aumenta la carga computacional
 */
void aumentar_carga() {
	double a=1.0, b=1.0;
	for (int i=0; i < 1000000; i++) a=a*b;
}

/**
 * Procedimiento que evalua la aptitud de un estudiante, 
 * calculando el valor correspondiente de la funcion objetivo
 * f(x)=x1*0,5+x2*-0,5+x3*0,75+x4*-0,75+x5*1
 * @param student: estudiante a evaluar
 */
void evaluate_student(struct Student* student) {
    float fx = student->solution[0] * 0.5
                    + student->solution[1] * -0.5
                    + student->solution[2] * 0.75
                    + student->solution[3] * -0.75
                    + student->solution[4] * 1.0;
    student->fx = fx;
    aumentar_carga();
    #pragma omp atomic
    num_evaluations++;   
     
}

void calculate_fitness_foreach_student(Student* classroom, int max_students){
     // Recorremos todos los estudiantes
    #pragma omp for
    for (int i = 0; i < max_students; i++ ) {
        evaluate_student(&classroom[i]);
    }
}

void init_student(Student *student) {
    // Recorrer las variables de decisión
    // unsigned int randomState = clock();
    for (int i = 0; i < MAX_SUBJECT; i++ )
        // student->solution[i]=rand_r(&randomState) % (MAX_NUMBER_SOLUTION + 1);
        student->solution[i]= dis(gen);
}


void print_classroom(Student classroom[], int max_students){
    // Recorremos todos los estudiantes
    for (int i = 0; i < max_students; i++ ) {
        print_student(&classroom[i], i);
        cout << endl;
    }
}

void init_classroom(Student* classroom, int max_students){
    // Recorremos todos los estudiantes
    #pragma omp for
    for (int i = 0; i < max_students; i++ ) {
        init_student(&classroom[i]);
    }
}

// Generar un número aleatorio en el rango [min, max]
int random_range(int min, int max) {
    unsigned int randomState = clock();
    return min + dis(gen);
}

// Algoritmo TLBO
void tlbo(Student classroom[], int max_students) {

    // para mejorar
    // ver generacion random, poner random especializado para paralizacion
    // ver la generacion de la primera poblacion que sea igual sin utilizar dos regiones paralizadas

    Student* teacher, newStudent;
    teacher = &classroom[0];

    int TF;
    int gen;

    #pragma omp parallel private(gen)
    {
        // Step3: For Iter in range(max_iter):  # loop max_iter times  
        for (gen = 0; gen < NUM_GENERATIONS && num_evaluations < MAX_EVALUATIONS; gen++){

            // ----------------------------------------------- TEACHING PHASE (teacher learn student)

            // calculate the mean of the classroom (mean of the each subject)
            int Xmean[MAX_SUBJECT];
            int student;

            #pragma omp for private(student)
            for (int subject = 0; subject < MAX_SUBJECT; subject++) {
                Xmean[subject] = 0;
                for (student = 0; student < max_students; student++){
                    Xmean[subject] += classroom[student].solution[subject]; // calculate mean
                    if(classroom[student].fx > teacher->fx)
                        teacher = &classroom[student]; // select the teacher (student with the best fitness)

                }
                Xmean[subject] /= max_students;
            }

            #pragma omp single
            // calculate teacher factor (TF)
            TF = random_range(1, 2);

            // calculate new solution in each student with the teacher lesson
            #pragma omp for
            for (int student = 0; student < max_students; student++){
                
                Student* newStudent = new Student();

                // calculate new solution
                int subject;
                for (subject = 0; subject < MAX_SUBJECT; subject++){
                    int r = random_range(0, 1);    
                    newStudent->solution[subject] = classroom[student].solution[subject] + r*(teacher->solution[subject] - TF*Xmean[subject]);
                    
                    // normalize solution
                    if(newStudent->solution[subject] > MAX_NUMBER_SOLUTION) newStudent->solution[subject] = MAX_NUMBER_SOLUTION;
                    if(newStudent->solution[subject] < MIN_NUMBER_SOLUTION) newStudent->solution[subject] = MIN_NUMBER_SOLUTION;
                }

                // evaluate new solution
                evaluate_student(newStudent);

                // #pragma omp critical
                // if (newStudent->fx == 20.25 && num_evaluations_found_solution == 0) {
                //     num_evaluations_found_solution = num_evaluations;
                // }

                // if new solution is better than the older, save the new solution
                if(newStudent->fx > classroom[student].fx) {
                    classroom[student] = *newStudent;

                    // ----------------------------------------------- LEARNING PHASE (student learn student)

                    // select another student to learn from
                    int another_student1 = random_range(0, max_students - 1);
                    while (another_student1 == student) 
                        another_student1 = random_range(0, max_students - 1);
                    Student* classmate1 = &classroom[another_student1];

                    int another_student2 = random_range(0, max_students - 1);
                    while (another_student2 == student) 
                        another_student2 = random_range(0, max_students - 1);
                    Student* classmate2 = &classroom[another_student2];

                    newStudent = new Student();

                    // calculate new solution
                    for (int subject = 0; subject < MAX_SUBJECT; subject++){
                        int r = random_range(0, 1);   

                        if(classmate1->solution[subject] > classmate2->solution[subject])
                            newStudent->solution[subject] = classroom[student].solution[subject] + r*(classmate1->solution[subject] - classmate2->solution[subject]);
                        else
                            newStudent->solution[subject] = classroom[student].solution[subject] + r*(classmate1->solution[subject] - classmate2->solution[subject]);

                        // normalize solution
                        if(newStudent->solution[subject] > MAX_NUMBER_SOLUTION) newStudent->solution[subject] = MAX_NUMBER_SOLUTION;
                        if(newStudent->solution[subject] < MIN_NUMBER_SOLUTION) newStudent->solution[subject] = MIN_NUMBER_SOLUTION;
                    }

                    evaluate_student(newStudent);
                    // #pragma omp critical
                    // if (newStudent->fx == 20.25 && num_evaluations_found_solution == 0) {
                    //     num_evaluations_found_solution = num_evaluations;
                    // }
                    
                    // if new solution is better than the older, save the new solution
                    if(newStudent->fx > classroom[student].fx){
                        classroom[student] = *newStudent;
                    }
                }
            }
            
        }
    }

    cout << "end with poblation with " << max_students << " students"<< endl;
    print_classroom(classroom, max_students);
    cout << "end with teacher " << endl;
    print_student(teacher, 0);
    cout << endl;
    cout << "---------------------------------------------------------" << endl;

}

double execWithThreads(Student classroom[], int max_students, int numThreads){
    omp_set_num_threads(numThreads);
    double Tini, Tfin = 0;
    num_evaluations = 0; // poner el contador 0 en cada ejecución
    num_evaluations_found_solution = 0; // poner el contador 0 en cada ejecución

    // copy classroom
    Student* copyClassroom = new Student[max_students];
    for (int i = 0; i < max_students; i++) copyClassroom[i] = classroom[i];

    Tini = omp_get_wtime(); // Obtener el tiempo inicial    
    cout << "init with poblation with " << max_students << " students"<< endl;
    print_classroom(copyClassroom, max_students);
    tlbo(copyClassroom, max_students); // Calcular el programa de los estudiantes
    Tfin = omp_get_wtime();  // Obtener el tiempo final y acumular el tiempo tardado
    cout << "TIME WITH " << numThreads << " THREAD: " << 1000 * (Tfin-Tini) << " milliseconds with " <<  num_evaluations << " evaluations"<< endl;
    cout << "SOLUTION FOUND IN " << num_evaluations_found_solution << " EVALUATIONS" << endl;
    return Tfin-Tini;
}

int main() {
    srand(time(NULL));

    // Declare classroom array in main and initialize
    Student classroom[MAX_STUDENT];

    #pragma omp parallel
    {

    // Step1: Randomly initialize Class of N students Xi ( i=1, 2, …, n)
    init_classroom(classroom, MAX_STUDENT);

    // Step2: Compute fitness value of all the students
    calculate_fitness_foreach_student(classroom, MAX_STUDENT);

    }
    
    double time1Thread = execWithThreads(classroom, MAX_STUDENT, 1);

    int allThreads = 4;
    double timeAllThreads = execWithThreads(classroom, MAX_STUDENT, allThreads);

    cout << "SPEEDUP x" << time1Thread/timeAllThreads << endl;
    
    return 0;
}
