# Xustificacion das melloras
No documento inicial preséntanse unha serie de melloras a seren probadas. A continuación faráse unha revisión de ditas melloras e explicaráse por que foron ou non probadas.

* A solución inicial
* O valor inicial do parámetro de control (T0).
* O mecanismo de arrefriamento.
* A velocidade de arrefriamento.
* A xeración da veciñanza.
* O criterio de parada

## Solución inicial
A solución inicial foi modificada probando a inserción dunha función que crease unha solución inicial aleatoria. Como era de agardar, esta opción non mellorou a solución atopada senón que a empeorou (_non en gran medida_). O código de inicio de dita solución inicial é o seguinte:
```C
void generateInitialSolution(){
int index=1,actualValue;
  actualSolution[0]=0;
  actualSolution[DIMENSION]=0;
  printActualSolution();
  for(;index<DIMENSION;){
    actualValue =(int) (1 + fmod(floor(rand()*numberOfCities-1),numberOfCities));
    for(;isvalueinarray(actualValue,actualSolution,numberOfCities);){
      actualValue++;
      actualValue=(int) fmod(actualValue,numberOfCities);
    }
    actualSolution[index]=actualValue;
    index++;
  }
  printActualSolution();
  calculateInitialDistance();
  minimalDistance=actualDistance;
  tDistance=minimalDistance;
  TEMPERATURE_ZERO= (MU/(-log(PHI))*actualDistance);
  TEMPERATURE=TEMPERATURE_ZERO;
}
```

## O valor inicial do parámetro de control (T0).
O valor inicial do parámetro de control foi modificado, tanto mediante a modificación das variables  **φ** e  **µ** como a través da asignación directa dun número enteiro arbitrario. Foron probados á súa vez número da orde de 0.0000x ata números enteiros e, en ningún dos casos se apreciou unha mellora con respecto á opción inicial polo que foi esta a que se mantivo.

## O mecanismo de arrefriamento.

O mecanismo de arrefriamento si foi modificado pois, ata o de agora o arrefriamento tan só modificaba a velocidade de arrefriamento. O que se fixo foi engadir un reinicio voraz cando se realizase un arrefriamento de forma que se lle dese un novo impulso á solución. Este reinicio voraz atópase condicionado por unha matriz que penaliza as maiores ocurrencias de pares de intercambio desta forma as solucións que se nos repiten de forma constante e que provocan que saia a mesma conclusión unha e outra vez vense evitadas en gran medida.

## A velocidade de arrefriamento
A velocidade de arrefriamento foi reducida e daptada para que se axustase o máximo posible ao número de iteracións previstas para o algoritmo.
```C
TEMPERATURE=TEMPERATURE_ZERO/(1+((TEMPERATURE_ZERO-TEMPERATURE)/(MAXITERATIONS*TEMPERATURE_ZERO*TEMPERATURE))*cooled);
```

## A xeración da veciñanza.
A xeración da veciñanza non foi modificada posto que agora mesmo se atopan xerados todos e cada un dos veciños posibles e o retardo da aplicación é asumible polo que non é ñun movemento intelixente neste caso penalizar as posibilidades de atopar unha mellor solución a cambio de gañar un tempo que apenas sería apreciable.

## O criterio de parada
O criterio de parada foi modificado ata poñer unha realación entre o número de iteracións que ocorrían sen mellora algunha pero, tras diversas probas comprobouse que este método non era o máis adecuado xa que si ocorrían moitas veces as mesmas repeticións antes de que se producise un resultado correcto polo que se desestimou o cambio  e se mantivo a versión orixinal con iteracións. Tamén se probou a cambiar o número máximo de iteracións a 100000 pero o resultado tampouco mellorou.
