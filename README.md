# secret-sharing
TFG


## Project Sctructure
```
practica/                   # Carpeta del proyecto
│
├── controladores/          # Carpeta con los controladores para ambas entidades
│   ├── cita.php            # Controlador de la entidad CITA
│   └── usuario.php         # Controlador de la entidad USUARIO
│
├── css/                    # Carpeta con el código CSS
│   └── style.css           # CSS propio
│
├── inc/                    # Carpeta con los sripts de PHP 
│   ├── alertas.php         # Generación de mensajes y alertas
│   ├── finalizacion.php    # Fin de la configuración (BBDD, buffers) 
│   ├── funciones.php       # Funciones de utilidad (parseo, conversión hora-hueco)
│   ├── inicializacion.php  # Inicio de la configuración (BBDD, constantes, sesión, ...)
│   └── plantilla.php       # Funciones para crear la  cabecera y el pie de página
│
├── js/                     # Carpeta con el código JavaScript
│   └── script.js           # Script para añadir validaciones en el cliente y AJAX en las tablas de las citas
│
├── modelos/                # Carpeta con los modelos de ambas entidades
│   ├── cita.php            # Modela de la entidad CITA
│   └── usuario.php         # Modelo de la entidad USUARIO
│
├── sql/                    # SQL del proyecto
│   └── practica.sql        # Script para crear la base de datos
│
├── vistas/                 # Carpeta con las vistas de ambas entidades
│   ├── index.php           # Página que permite ver la agenda a los usuarios no identificados
│   ├── misCitas.php        # Página que permite al usuario registrado ver, filtrar y anular sus citas
│   ├── nuevaCita.php       # Página que permite al usuario registrado pedir una cita
│   ├── registro.php        # Pagina que contiene el formulario de registro en el sistema
│   ├── sesion.php          # Página que contiene el formulario de identificación en el sistema
│   └── usuario.php         # Página principal del usuario registrado
│
└── index.php               # Página principal del proyecto
```