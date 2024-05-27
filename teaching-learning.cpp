#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <sys/time.h>

#define MAX_SUBJECT 5
#define MAX_NUMBER_SOLUTION 9
#define MIN_NUMBER_SOLUTION 0
#define MAX_STUDENT 15
#define NUM_GENERATIONS 100

int num_evaluations = 0;

using namespace std;

//type student
struct Student {
    int solution[MAX_SUBJECT];
    float fx;
};

//Tipo de datos necesario para el qsort
typedef int (*compfn) (const void*, const void*);


//clase de estudiantes
Student classroom[MAX_STUDENT];

// Imprimir una solución con su índice, valores y aptitud
void print_student(Student *student, int index) {
    cout << "Student " << index << ": ";
    for (int i = 0; i < MAX_SUBJECT; ++i) {
        cout << student->solution[i] << " ";
    }
    cout << "| fx: " << student->fx;
}

/**
 * Funcion utilizada en el qsort para ordenar el enjambre en funcion del valor de fx
 * @param student1 
 * @param student2
 * @return 1, -1 o 0
 */
int compareFX(struct Student* student1, struct Student* student2) {
    if (student1->fx < student2->fx) return 1;
    if (student1->fx > student2->fx) return -1;
    return 0;
}

Student* select_the_teacher(){
    qsort(classroom, MAX_STUDENT, sizeof(Student), (compfn) compareFX);
    return &classroom[0]; // the first student is the best in the classroom
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
    num_evaluations++;
}

void calculate_fitness_foreach_student(){
     // Recorremos todos los estudiantes
    for (int i = 0; i < MAX_STUDENT; i++ ) {
        evaluate_student(&classroom[i]);
    }
}

void init_student(Student *student) {
    // Recorrer las variables de decisión
    for (int i = 0; i < MAX_SUBJECT; i++ )
        student->solution[i]=rand() % (MAX_NUMBER_SOLUTION + 1);
}

void init_classroom(){
    // Recorremos todos los estudiantes
    for (int i = 0; i < MAX_STUDENT; i++ ) {
        init_student(&classroom[i]);
    }
}


// Generar un número aleatorio en el rango [min, max]
int random_range(int min, int max) {
    return min + rand() % (max - min + 1);
}

// Algoritmo TLBO
void tlbo() {

    Student* teacher;

    // Step1: Randomly initialize Class of N students Xi ( i=1, 2, …, n)
    init_classroom();

    // Step2: Compute fitness value of all the students
    calculate_fitness_foreach_student();

    // Step3: For Iter in range(max_iter):  # loop max_iter times  
    for (int gen = 0; gen < NUM_GENERATIONS; gen++){

        cout << "START GEN " << gen + 1 << endl;

        // ----------------------------------------------- TEACHING PHASE (teacher learn student)

        // calculate the mean of the classroom (mean of the each subject)
        int Xmean[MAX_SUBJECT];

        for (int subject = 0; subject < MAX_SUBJECT; subject++) {
            Xmean[subject] = 0;
            for (int student = 0; student < MAX_STUDENT; student++) Xmean[subject] += classroom[student].solution[subject];
            Xmean[subject] /= MAX_STUDENT;
        }

        // print the mean of the class in each subject
        cout << "Xmean " << ": ";
        for (int i = 0; i < MAX_SUBJECT; ++i)
            cout << Xmean[i] << " ";
        cout << endl;


        // select the teacher (student with the best fitness)
        teacher = select_the_teacher();
        cout << "Teacher selected:" << ' ';
        print_student(teacher, 0);
        cout << endl;

        

        // calculate teacher factor (TF)
        int TF = random_range(1, 2);

        // calculate new solution in each student with the teacher leason
        for (int student = 0; student < MAX_STUDENT; student++){
            print_student(&classroom[student], student);
            cout << endl;
   
            // start the new student
            Student* newStudent = new Student();
            // calculate new solution
            for (int subject = 0; subject < MAX_SUBJECT; subject++){
                

                int r = random_range(0, 1);    
                newStudent->solution[subject] = classroom[student].solution[subject] + r*(teacher->solution[subject] - TF*Xmean[subject]);
                
                // normalice solution
                if(newStudent->solution[subject] > MAX_NUMBER_SOLUTION) newStudent->solution[subject] = MAX_NUMBER_SOLUTION;
                if(newStudent->solution[subject] < MIN_NUMBER_SOLUTION) newStudent->solution[subject] = MIN_NUMBER_SOLUTION;
            }

            // evalutate new solution
            evaluate_student(newStudent);
            // if new solution is better than the older, save the new solution
            if(newStudent->fx > classroom[student].fx) {
                
                cout << endl << "CHANGE STUDENT[" << student << "] TEACHING PHASE ";
                cout << "BEFORE [";
                print_student(&classroom[student], student);
                cout << "] -> AFTER[";
                print_student(newStudent, student);
                cout << ']' << endl;
                classroom[student] = *newStudent;

                // ----------------------------------------------- LEARNING PHASE (student learn student)

                // select anothers student to learn
                int another_student1 = random_range(0, MAX_STUDENT - 1);
                while (another_student1 == student) 
                    another_student1 = random_range(0, MAX_STUDENT - 1);
                Student* classmate1 = &classroom[another_student1];

                int another_student2 = random_range(0, MAX_STUDENT - 1);
                while (another_student2 == student) 
                    another_student2 = random_range(0, MAX_STUDENT - 1);
                Student* classmate2 = &classroom[another_student2];

                newStudent = new Student();

                // calculate new solution
                for (int subject = 0; subject < MAX_SUBJECT; subject++){
                    int r = random_range(0, 1);   

                    if(classmate1->solution[subject] > classmate2->solution[subject])
                        newStudent->solution[subject] = classroom[student].solution[subject] + 
                                                        r*(classmate1->solution[subject] - classmate2->solution[subject]);
                    else
                        newStudent->solution[subject] = classroom[student].solution[subject] 
                                                        + r*(classmate1->solution[subject] - classmate2->solution[subject]);

                    // normalice solution
                    if(newStudent->solution[subject] > MAX_NUMBER_SOLUTION) newStudent->solution[subject] = MAX_NUMBER_SOLUTION;
                    if(newStudent->solution[subject] < MIN_NUMBER_SOLUTION) newStudent->solution[subject] = MIN_NUMBER_SOLUTION;

                }

                evaluate_student(newStudent);
                // if new solution is better than the older, save the new solution
                if(newStudent->fx > classroom[student].fx){
                    cout << "CHANGE STUDENT[" << student << "] LEARNING PHASE ";
                    cout << "BEFORE [";
                    print_student(&classroom[student], student);
                    cout << "] -> AFTER[";
                    print_student(newStudent, student);
                    cout << ']' << endl;
                    classroom[student] = *newStudent;
                }
            }
           

        }

    // Step 4: Return best student from class
    cout << "BEST SOLUTION IN GEN " << gen + 1<< ' ';
    print_student(teacher, 0);
    cout << endl;
    cout << "------------------------------------------" << endl;
    if (teacher->fx == 20.25) break;


    }

// Step 4: Return best student from class
cout << "BEST SOLUTION FOUND" << ' ';
print_student(teacher, 0);
cout << endl;
cout << "------------------------------------------" << endl;

cout << "NUM EVALUATIONS : " << num_evaluations <<  endl;

}

int main() {
    srand(time(NULL));
    tlbo(); 
    return 0;
}