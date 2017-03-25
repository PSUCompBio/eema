CC = g++
CFLAGS = -g -Wall
INCLUDES = -I./third-party-libs/eigen3 -I./headers

# If you increase the number of source subdirectories, then add include the new name here and
# add the new rule at the end.
# or use the below option i.e., ./source/*.
SRC_0 = ./source
SRC_1 = ./source/*
SRC_2 = ./source/*/*
SRC_3 = ./source/*/*/*
HDR = ./headers
OBJ = ./obj
MAIN = main
JOB = ./examples/example-1

OBJECTS := $(OBJ)/$(MAIN).o $(OBJ)/fe_vector2text.o\
			 $(OBJ)/fe_apply_bc.o\
			 $(OBJ)/fe_text2matrix.o\
			 $(OBJ)/fe_matrix2text.o\
			 $(OBJ)/fe_assemble_mass.o\
			 $(OBJ)/fe_calArea_4.o\
			 $(OBJ)/fe_calVolume.o\
			 $(OBJ)/fe_calJacobian.o\
			 $(OBJ)/fe_calSimpTransformation.o\
			 $(OBJ)/fe_calTransformation.o\
			 $(OBJ)/fe_concatenate_vector2matrix.o\
			 $(OBJ)/fe_display.o\
			 $(OBJ)/fe_dn_actual_8.o\
			 $(OBJ)/fe_dn_iso_8.o\
			 $(OBJ)/fe_get_mats.o\
			 $(OBJ)/fe_getforce.o\
			 $(OBJ)/fe_guass_points.o\
			 $(OBJ)/fe_guass_points_3d.o\
			 $(OBJ)/fe_guass_weights.o\
			 $(OBJ)/fe_guass_weights_3d.o\
			 $(OBJ)/fe_mass_hex.o\
			 $(OBJ)/fe_shapeMatrix.o\
			 $(OBJ)/fe_shapes.o\
			 $(OBJ)/fe_strDispMatrix.o\
			 $(OBJ)/fe_stressUpdate.o\
			 $(OBJ)/fe_stressUpdate_1d.o\
			 $(OBJ)/fe_getTimeStep.o\
			 $(OBJ)/fe_calTimeStep.o\
			 $(OBJ)/fe_calWaveSpeed.o\
			 $(OBJ)/fe_transformMass.o\
			 $(OBJ)/fe_calDefGrad.o\
			 $(OBJ)/fe_write.o\
			 $(OBJ)/Mesh.o\
			 $(OBJ)/fe_vtk.o\
			 $(OBJ)/fe_vtu.o\
			 $(OBJ)/fe_update.o\
			 $(OBJ)/Materials.o\
			 $(OBJ)/fe_mainEXPLICIT.o\
			 $(OBJ)/fe_find_index.o\
			 $(OBJ)/fe_gather.o\
			 $(OBJ)/fe_scatter.o\
			 $(OBJ)/fe_tensor2voigt.o\
			 $(OBJ)/fe_voigt2tensor.o\
			 $(OBJ)/fe_find.o\
			 $(OBJ)/fe_mainRead.o\
			 $(OBJ)/fe_function.o\
			 $(OBJ)/fe_mooneyrivlin_hyperelastic.o\
			 $(OBJ)/fe_ogden_hyperelastic.o\
			 $(OBJ)/fe_simple_elastic.o\
			 $(OBJ)/fe_saintvenant_elastic.o\

all: $(MAIN)

$(MAIN):	$(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJECTS) -o $(MAIN)

$(OBJ)/$(MAIN).o:	$(MAIN).cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $(MAIN).cpp -o $(OBJ)/$(MAIN).o

$(OBJ)/%.o: $(SRC_0)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ)/%.o: $(SRC_1)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ)/%.o: $(SRC_2)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ)/%.o: $(SRC_3)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(MAIN) obj/*.o *~

res:
	./clean.sh
