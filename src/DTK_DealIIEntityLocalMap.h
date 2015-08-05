#ifndef DTK_DEALIIENTITYLOCALMAP_HPP
#define DTK_DEALIIENTITYLOCALMAP_HPP

#include <DTK_EntityLocalMap.hpp>

#include <deal.II/fe/mapping.h>

template <int dim,int spacedim>
class DealIIEntityLocalMap : public DataTransferKit::EntityLocalMap
{
  public:

    DealIIEntityLocalMap(std::shared_ptr<dealii::Mapping<dim,spacedim> const> mapping);

    /*
     * \brief Set parameters for mapping.
     * \param parameters Parameters for mapping.
     */
    void setParameters( const Teuchos::ParameterList& parameters ) override;

    /*!
     * \brief Return the entity measure with respect to the parameteric
     * dimension (volume for a 3D entity, area for 2D, and length for 1D).
     * \param entity Compute the measure for this entity.
     * \return The measure of the entity.
     */
    double measure( const DataTransferKit::Entity& entity ) const override;

    /*!
     * \brief Return the centroid of the entity.
     * \param centroid A view of the centroid coordinates. This view will
     * be allocated. Assign a view of your centroid to this view.
     */
    void centroid( const DataTransferKit::Entity& entity,
       const Teuchos::ArrayView<double>& centroid ) const override;

    /*!
     * \brief (Safeguard the reverse map) Perform a safeguard check for
     * mapping a point to the reference space of an entity using the given
     * tolerance.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the safeguard check.
     * \param physical_point A view into an array of size physicalDimension()
     * containing the coordinates of the point to map.
     * \return Return true if it is safe to map to the reference frame.
     */
    bool isSafeToMapToReferenceFrame(
  const DataTransferKit::Entity& entity,
  const Teuchos::ArrayView<const double>& physical_point ) const override;

    /*!
     * \brief (Reverse Map) Map a point to the reference space of an
     * entity. Return the parameterized point.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the mapping procedure.
     * \param physical_point A view into an array of size physicalDimension()
     * containing the coordinates of the point to map.
     * \param reference_point A view into an array of size physicalDimension()
     * to write the reference coordinates of the mapped point.
     * \return Return true if the map to reference frame succeeded.
     */
    bool mapToReferenceFrame( 
        const DataTransferKit::Entity& entity,
        const Teuchos::ArrayView<const double>& physical_point,
        const Teuchos::ArrayView<double>& reference_point ) const override;

    /*!  
     * \brief Determine if a reference point is in the parameterized space of
     * an entity.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the point inclusion check.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \return True if the point is in the reference space, false if not.
     */
    bool checkPointInclusion( 
        const DataTransferKit::Entity& entity,
        const Teuchos::ArrayView<const double>& reference_point ) const override;

    /*!
     * \brief (Forward Map) Map a reference point to the physical space of an
     * entity. 
     * \param entity Perfrom the mapping for this entity.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \param physical_point A view into an array of size physicalDimension()
     * to write the coordinates of physical point.
     */
    void mapToPhysicalFrame( 
        const DataTransferKit::Entity& entity,
        const Teuchos::ArrayView<const double>& reference_point,
        const Teuchos::ArrayView<double>& physical_point ) const override;

    /*!
     * \brief Compute the normal on a face (3D) or edge (2D) at a given
     * reference point. A default implementation is provided using a finite
     * difference scheme.
     * \param entity Compute the normal for this entity.
     * \param parent_entity The adjacent parent entity used to determine which
     * direction is outward. The parent entity should be of a higher
     * topological dimension than the entity and be adjacent to the entity.
     * \param reference_point Compute the normal at this reference point.
     * \param normal A view into an array of size physicalDimension() to write
     * the normal.
     */
    void normalAtReferencePoint( 
        const DataTransferKit::Entity& entity,
        const DataTransferKit::Entity& parent_entity,
        const Teuchos::ArrayView<const double>& reference_point,
        const Teuchos::ArrayView<double>& normal ) const override;

  private:
    std::shared_ptr<dealii::Mapping<dim,spacedim> const> dealii_mapping;
};

#endif
