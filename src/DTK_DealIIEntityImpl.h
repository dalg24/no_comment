#ifndef DTK_DEALIIENTITYIMPL_HPP
#define DTK_DEALIIENTITYIMPL_HPP

#include <no_comment/DTK_DealIIEntityExtraData.h>
#include <DTK_EntityImpl.hpp>
#include <DTK_Types.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

template <int structdim,int dim,int spacedim>
class DealIIEntityImpl : public DataTransferKit::EntityImpl
{
public:

    DealIIEntityImpl(dealii::TriaAccessor<structdim,dim,spacedim> const & tria_accessor,
        Teuchos::RCP<std::vector<std::set<
        typename dealii::Triangulation<dim,spacedim>::active_cell_iterator>>> vertex_to_cell,
        Teuchos::RCP<std::map<unsigned int, unsigned long long int>> local_to_global_vertex_id);

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    DataTransferKit::EntityId id() const override;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const override;

    /*!
     * \brief Return the topological dimension of the entity.  
     *
     * \return The topological dimension of the entity. Any parametric
     * coordinates describing the entity will be of this dimension.
     */
    int topologicalDimension() const override;

    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    int physicalDimension() const override;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const override;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const override;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const override;

    /*!
     * \brief Get the extra data on the entity.
     */
    Teuchos::RCP<DataTransferKit::EntityExtraData> extraData() const override;

    /*!
     * \brief Provide a one line description of the object.
     */
    std::string description() const override
    { return std::string("deal.II Entity"); }

    /*!
     * \brief Provide a verbose description of the object.
     */
    void describe(
      	Teuchos::FancyOStream& out,
      	const Teuchos::EVerbosityLevel verb_level ) const override;

private:
    Teuchos::RCP<DealIIEntityExtraData<structdim,dim,spacedim>> extra_data;
};

#endif
