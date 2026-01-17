/*
Autor: Ivan Rugerio
Titulo: Sistema de Gestión de Préstamos y Análisis de Riesgo Geoespacial (PostGIS)
Descripción: Base de datos para administración de solicitudes bancarias con análisis espacial de riesgos en CDMX.
Tecnologías: PostgreSQL, PostGIS
Version: 1.0
*/

-- --- 1. CONFIGURACIÓN DE LA BASE DE DATOS ---

-- CREAR BD (Ejecutar solo si no existe)
-- CREATE DATABASE ir_pf
-- 	WITH
-- 	OWNER = postgres
-- 	ENCODING = 'UTF8';

-- CONECTAR A BD
-- \c ir_pf

-- --- 2. CREACIÓN DE TABLAS ---

CREATE TABLE personas (
	id_curp VARCHAR(18) NOT NULL,
	nombre VARCHAR(100), 
	apellido_paterno VARCHAR(100),
	apellido_materno VARCHAR(100),
	fecha_nacimiento DATE,
	edad INTEGER,
	correo_electronico VARCHAR(100),
	telefono BIGINT,
	delegacion_id INTEGER,
	escuela_id INTEGER,
	CONSTRAINT pk_curp PRIMARY KEY(id_curp)
);

CREATE TABLE delegacion (
	id_delegacion INTEGER NOT NULL,
	nombre_delegacion VARCHAR(100), 
	CONSTRAINT pk_delegacion PRIMARY KEY(id_delegacion)
);

CREATE TABLE escuela (
	id_escuela INTEGER NOT NULL,
	nombre_escuela VARCHAR(100), 
	ubicacion VARCHAR(100),
	CONSTRAINT pk_escuela PRIMARY KEY(id_escuela)
);

CREATE TABLE prestamo (
	id_prestamo INTEGER NOT NULL,
	fecha_inicio DATE,
	fecha_final DATE,
	monto INTEGER,
	cantidad_prestamos INTEGER,
	motivo VARCHAR(100),
	curp_id VARCHAR(18),
	banco_id INTEGER,
	CONSTRAINT pk_prestamo PRIMARY KEY(id_prestamo)
);

CREATE TABLE banco (
	id_banco INTEGER NOT NULL,
	nombre_banco VARCHAR(100), 
	CONSTRAINT pk_banco PRIMARY KEY(id_banco)
);

CREATE TABLE solicitud (
	id_solicitud INTEGER NOT NULL,
	fecha_solicitud DATE,
	factible CHAR(2),
	banco_id INTEGER,
	curp_id VARCHAR(18),
	latitud DOUBLE PRECISION,
	longitud DOUBLE PRECISION,
	CONSTRAINT pk_solicitud PRIMARY KEY(id_solicitud)
);

-- --- 3. RELACIONES (LLAVES FORÁNEAS) ---

ALTER TABLE personas ADD CONSTRAINT fk_delegacion FOREIGN KEY (delegacion_id)
REFERENCES delegacion (id_delegacion) ON DELETE RESTRICT ON UPDATE CASCADE;

ALTER TABLE personas ADD CONSTRAINT fk_escuela FOREIGN KEY (escuela_id)
REFERENCES escuela (id_escuela) ON DELETE RESTRICT ON UPDATE CASCADE;  

ALTER TABLE prestamo ADD CONSTRAINT fk_curp_prestamo FOREIGN KEY (curp_id)
REFERENCES personas (id_curp) ON DELETE RESTRICT ON UPDATE CASCADE;  

ALTER TABLE prestamo ADD CONSTRAINT fk_banco_prestamo FOREIGN KEY (banco_id)
REFERENCES banco (id_banco) ON DELETE RESTRICT ON UPDATE CASCADE;  

ALTER TABLE solicitud ADD CONSTRAINT fk_curp_solicitud FOREIGN KEY (curp_id)
REFERENCES personas (id_curp) ON DELETE RESTRICT ON UPDATE CASCADE;  

ALTER TABLE solicitud ADD CONSTRAINT fk_banco_solicitud FOREIGN KEY (banco_id)
REFERENCES banco (id_banco) ON DELETE RESTRICT ON UPDATE CASCADE;  

-- --- 4. POBLADO DE DATOS (CATÁLOGOS) ---

INSERT INTO delegacion(id_delegacion,nombre_delegacion) VALUES
	(1,'Azcapotzalco'), (2,'Coyoacan'), (3,'Cuajimalpa de Morelos'),
	(4,'Gustavo A. Madero'), (5,'Iztacalco'), (6,'Iztapalapa'),
    (7,'La Magdalena Contreras'), (8,'Milpa Alta'), (9,'Alvaro Obregon'),
	(10,'Tlahuac'), (11,'Tlalpan'), (12,'Xochimilco'),
    (13,'Benito Juarez'), (14,'Cuauhtemoc'), (15,'Miguel Hidalgo'), 
	(16,'Venustiano Carranza');

INSERT INTO escuela(id_escuela,nombre_escuela,ubicacion) VALUES
	(1,'Facultad de Arquitectura','Circuito Escolar s/n, CU'),
	(2,'Facultad de Artes y Diseno','Xochimilco'),
	(3,'Facultad de Ciencias','Circuito Exterior s/n, CU'),
	(4,'Facultad de Ciencias Politicas y Sociales','Circuito Mario de la Cueva s/n, CU'),
	(5,'Facultad de Contaduria y Administracion','Av. Universidad 3000, CU'),
	(10,'Facultad de Ingenieria','Circuito Escolar s/n, CU'),
	(11,'Facultad de Medicina','Circuito Interior s/n, CU'),
    -- (Se simplificó la lista por brevedad, agrega el resto si es necesario)
	(23,'Escuela Nacional de Ciencias Forenses','Av. Universidad 3000, CU');

INSERT INTO banco(id_banco,nombre_banco) VALUES
	(1,'BBVA'), (2,'Santander'), (3,'Banamex'), (4,'Banco Azteca'), (5,'Nu');

-- --- 5. CARGA DE DATOS MASIVOS (CSV) ---
-- NOTA: Estas líneas están comentadas para portabilidad. 
-- Descomentar y ajustar la ruta si se tienen los archivos CSV localmente.

-- \copy personas FROM './datos/personas.csv' DELIMITER ',' CSV HEADER;
-- \copy solicitud FROM './datos/solicitudes.csv' DELIMITER ',' CSV HEADER;
-- \copy prestamo FROM './datos/prestamos.csv' DELIMITER ',' CSV HEADER;


-- --- 6. GEOMETRÍA Y ANÁLISIS ESPACIAL (POSTGIS) ---

CREATE EXTENSION IF NOT EXISTS postgis;

-- Añadir geometría a delegaciones (Se requiere tabla externa de polígonos)
ALTER TABLE delegacion ADD COLUMN geom geometry(MULTIPOLYGON, 4326);

-- UPDATE delegacion d SET geom = dt.geom
-- FROM poligonos_alcaldias_cdmx dt WHERE d.id_delegacion = dt.id; 
-- DROP TABLE poligonos_alcaldias_cdmx;

-- Generar puntos espaciales para las solicitudes
ALTER TABLE solicitud ADD COLUMN geom geometry(Point, 4326);

UPDATE solicitud
SET geom = ST_SetSRID(ST_MakePoint(longitud, latitud), 4326);

-- --- 7. ANÁLISIS DE RIESGOS (INCIDENCIAS) ---

CREATE TABLE incidencias (
    id SERIAL PRIMARY KEY,
    descripcion TEXT,
    geom GEOMETRY(POLYGON, 4326)
);

-- Generación de zonas de riesgo (Buffers de 1km)
INSERT INTO incidencias (descripcion, geom) VALUES 
    ('Incidencia 1 - Vallejo', ST_Transform(ST_Buffer(ST_Transform(ST_SetSRID(ST_MakePoint(-99.14304224, 19.45430286), 4326), 3857), 1000), 4326)),
    ('Incidencia 2 - Zona Centro', ST_Transform(ST_Buffer(ST_Transform(ST_SetSRID(ST_MakePoint(-99.10908567, 19.44545983), 4326), 3857), 1000), 4326)),
    ('Incidencia 3 - San Jeronimo', ST_Transform(ST_Buffer(ST_Transform(ST_SetSRID(ST_MakePoint(-99.22515584, 19.36542339), 4326), 3857), 1000), 4326)),
    ('Incidencia 4 - Tlalpan', ST_Transform(ST_Buffer(ST_Transform(ST_SetSRID(ST_MakePoint(-99.18513762, 19.28590667), 4326), 3857), 1000), 4326)),
    ('Incidencia 5 - Iztapalapa', ST_Transform(ST_Buffer(ST_Transform(ST_SetSRID(ST_MakePoint(-99.07651673, 19.30721507), 4326), 3857), 1000), 4326));

-- Identificación de solicitudes en zonas de riesgo
-- (Requiere tabla 'indi_pv' existente en la BD)
/*
CREATE TABLE intersecciones AS
SELECT ai.id AS area_id, ip.id AS indi_pv_id, ST_Intersection(ai.geom, ip.geom) AS geom
FROM incidencias ai JOIN indi_pv ip ON ST_Intersects(ai.geom, ip.geom);

CREATE TABLE solicitudes_en_zonas_riesgo AS
SELECT s.id_solicitud, s.curp_id, s.geom AS solicitud_geom, i.area_id, i.geom AS zona_interseccion
FROM solicitud s JOIN intersecciones i ON ST_Intersects(s.geom, i.geom);
*/

-- --- 8. EJEMPLO DE RUTAS ---

CREATE TABLE sitios_ruta (
    id SERIAL PRIMARY KEY,
    nombre TEXT,
    latitud DOUBLE PRECISION,
    longitud DOUBLE PRECISION,
    geom GEOMETRY(POINT, 4326)
);

INSERT INTO sitios_ruta (nombre, latitud, longitud, geom) VALUES
('CDMX', 19.434996, -99.141227, ST_SetSRID(ST_MakePoint(-99.141227, 19.434996), 4326)),
('Puerto Vallarta', 20.686059, -105.2919229, ST_SetSRID(ST_MakePoint(-105.2919229, 20.686059), 4326)),
('Cartagena', 10.3873232, -75.5401143, ST_SetSRID(ST_MakePoint(-75.5401143, 10.3873232), 4326));

CREATE TABLE ruta_sitios AS
SELECT 1 AS id, ST_MakeLine(geom ORDER BY id) AS geom FROM sitios_ruta;
